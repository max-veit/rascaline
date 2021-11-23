use std::os::raw::c_char;
use std::ffi::CStr;
use std::ops::{Deref, DerefMut};

use rascaline::{Calculator, System, Error, CalculationOptions, SelectedIndexes};
use rascaline::descriptor::IndexesBuilder;

use super::utils::copy_str_to_c;
use super::{catch_unwind, rascal_status_t};

use super::descriptor::{rascal_descriptor_t, rascal_indexes_t};
use super::system::rascal_system_t;

/// Opaque type representing a `Calculator`
#[allow(non_camel_case_types)]
pub struct rascal_calculator_t(Calculator);

impl Deref for rascal_calculator_t {
    type Target = Calculator;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for rascal_calculator_t {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/// Create a new calculator with the given `name` and `parameters`.
///
/// @verbatim embed:rst:leading-asterisk
///
/// The list of available calculators and the corresponding parameters are in
/// the :ref:`main documentation <calculators-list>`. The ``parameters`` should
/// be formatted as JSON, according to the requested calculator schema.
///
/// @endverbatim
///
/// All memory allocated by this function can be released using
/// `rascal_calculator_free`.
///
/// @param name name of the calculator as a NULL-terminated string
/// @param parameters hyper-parameters of the calculator, JSON-formatted in a
///                   NULL-terminated string
///
/// @returns A pointer to the newly allocated calculator, or a `NULL` pointer in
///          case of error. In case of error, you can use `rascal_last_error()`
///          to get the error message.
#[no_mangle]
#[allow(clippy::module_name_repetitions)]
pub unsafe extern fn rascal_calculator(name: *const c_char, parameters: *const c_char) -> *mut rascal_calculator_t {
    let mut raw = std::ptr::null_mut();
    let unwind_wrapper = std::panic::AssertUnwindSafe(&mut raw);
    let status = catch_unwind(move || {
        check_pointers!(name, parameters);
        let name = CStr::from_ptr(name).to_str()?;
        let parameters = CStr::from_ptr(parameters).to_str()?;
        let calculator = Calculator::new(name, parameters.to_owned())?;
        let boxed = Box::new(rascal_calculator_t(calculator));

        *unwind_wrapper.0 = Box::into_raw(boxed);
        Ok(())
    });

    if !status.is_success() {
        return std::ptr::null_mut();
    }

    return raw;
}

/// Free the memory associated with a `calculator` previously created with
/// `rascal_calculator`.
///
/// If `calculator` is `NULL`, this function does nothing.
///
/// @param calculator pointer to an existing calculator, or `NULL`
///
/// @returns The status code of this operation. If the status is not
///          `RASCAL_SUCCESS`, you can use `rascal_last_error()` to get the
///          full error message.
#[no_mangle]
pub unsafe extern fn rascal_calculator_free(calculator: *mut rascal_calculator_t) -> rascal_status_t {
    catch_unwind(|| {
        if !calculator.is_null() {
            let boxed = Box::from_raw(calculator);
            std::mem::drop(boxed);
        }

        Ok(())
    })
}

/// Get a copy of the name of this calculator in the `name` buffer of size
/// `bufflen`.
///
/// `name` will be NULL-terminated by this function. If the buffer is too small
/// to fit the whole name, this function will return
/// `RASCAL_BUFFER_SIZE_ERROR`
///
/// @param calculator pointer to an existing calculator
/// @param name string buffer to fill with the calculator name
/// @param bufflen number of characters available in the buffer
///
/// @returns The status code of this operation. If the status is not
///          `RASCAL_SUCCESS`, you can use `rascal_last_error()` to get the full
///          error message.
#[no_mangle]
pub unsafe extern fn rascal_calculator_name(
    calculator: *const rascal_calculator_t,
    name: *mut c_char,
    bufflen: usize
) -> rascal_status_t {
    catch_unwind(|| {
        check_pointers!(calculator, name);
        copy_str_to_c(&(*calculator).name(), name, bufflen)?;
        Ok(())
    })
}

/// Get a copy of the parameters used to create this calculator in the
/// `parameters` buffer of size `bufflen`.
///
/// `parameters` will be NULL-terminated by this function. If the buffer is too
/// small to fit the whole name, this function will return
/// `RASCAL_BUFFER_SIZE_ERROR`.
///
/// @param calculator pointer to an existing calculator
/// @param parameters string buffer to fill with the parameters used to create
///                   this calculator
/// @param bufflen number of characters available in the buffer
///
/// @returns The status code of this operation. If the status is not
///          `RASCAL_SUCCESS`, you can use `rascal_last_error()` to get the full
///          error message.
#[no_mangle]
pub unsafe extern fn rascal_calculator_parameters(
    calculator: *const rascal_calculator_t,
    parameters: *mut c_char,
    bufflen: usize
) -> rascal_status_t {
    catch_unwind(|| {
        check_pointers!(calculator, parameters);
        copy_str_to_c((*calculator).parameters(), parameters, bufflen)?;
        Ok(())
    })
}

/// Get the default number of features this `calculator` will produce in the
/// `features` parameter.
///
/// This number corresponds to the size of second dimension of the `values` and
/// `gradients` arrays in the `rascal_descriptor_t` after a call to
/// `rascal_calculator_compute`.
///
/// @param calculator pointer to an existing calculator
/// @param features pointer to an integer to be filled with the number of features
///
/// @returns The status code of this operation. If the status is not
///          `RASCAL_SUCCESS`, you can use `rascal_last_error()` to get the full
///          error message.
#[no_mangle]
pub unsafe extern fn rascal_calculator_features_count(
    calculator: *const rascal_calculator_t,
    features: *mut usize
) -> rascal_status_t {
    catch_unwind(|| {
        check_pointers!(calculator, features);
        *features = (*calculator).default_features().count();
        Ok(())
    })
}

/// Options that can be set to change how a calculator operates.
#[repr(C)]
pub struct rascal_calculation_options_t {
    /// Copy the data from systems into native `SimpleSystem`. This can be
    /// faster than having to cross the FFI boundary too often.
    use_native_system: bool,
    /// List of samples on which to run the calculation. You can set
    /// `selected_samples.values` to `NULL` to run the calculation on all
    /// samples. If necessary, gradients samples will be derived from the
    /// values given in selected_samples.
    selected_samples: rascal_indexes_t,
    /// List of features on which to run the calculation. You can set
    /// `selected_features.values` to `NULL` to run the calculation on all
    /// features.
    selected_features: rascal_indexes_t,
}

fn selected_indexes(selected: &rascal_indexes_t) -> Result<SelectedIndexes, Error> {
    if selected.values.is_null() {
        return Ok(SelectedIndexes::All);
    }
    let names = unsafe {
        std::slice::from_raw_parts(selected.names, selected.size)
    };
    let names = names.iter()
        .map(|&s| unsafe {
            CStr::from_ptr(s).to_str()
        })
        .collect::<Result<Vec<_>, _>>()?;

    let values = unsafe {
        std::slice::from_raw_parts(selected.values.cast(), selected.size * selected.count)
    };

    let mut builder = IndexesBuilder::new(names);
    for chunk in values.chunks(selected.size) {
        builder.add(chunk);
    }

    return Ok(SelectedIndexes::Some(builder.finish()));
}

#[allow(clippy::doc_markdown)]
/// Run a calculation with the given `calculator` on the given `systems`,
/// storing the resulting data in the `descriptor`.
///
/// @param calculator pointer to an existing calculator
/// @param descriptor pointer to an existing descriptor for data storage
/// @param systems pointer to an array of systems implementation
/// @param systems_count number of systems in `systems`
/// @param options options for this calculation
///
/// @returns The status code of this operation. If the status is not
///          `RASCAL_SUCCESS`, you can use `rascal_last_error()` to get the full
///          error message.
#[no_mangle]
pub unsafe extern fn rascal_calculator_compute(
    calculator: *mut rascal_calculator_t,
    descriptor: *mut rascal_descriptor_t,
    systems: *mut rascal_system_t,
    systems_count: usize,
    options: rascal_calculation_options_t,
) -> rascal_status_t {
    catch_unwind(|| {
        if systems_count == 0 {
            log::warn!("0 systems given to rascal_calculator_compute, we will do nothing");
            return Ok(());
        }
        check_pointers!(calculator, descriptor, systems);

        // Create a Vec<Box<dyn System>> from the passed systems
        let c_systems = std::slice::from_raw_parts_mut(systems, systems_count);
        let mut systems = Vec::with_capacity(c_systems.len());
        for system in c_systems {
            systems.push(Box::new(system) as Box<dyn System>);
        }

        let options = CalculationOptions {
            use_native_system: options.use_native_system,
            selected_samples: selected_indexes(&options.selected_samples)?,
            selected_features: selected_indexes(&options.selected_features)?,
        };

        (*calculator).compute(&mut systems, &mut *descriptor, options)
    })
}
