use std::collections::BTreeSet;

use ndarray::parallel::prelude::*;

use crate::descriptor::{SamplesBuilder, IndexValue, Indexes, IndexesBuilder};
use crate::descriptor::TwoBodiesSpeciesSamples;

use crate::{CalculationOptions, Calculator, SelectedIndexes};
use crate::{Descriptor, Error, System};

use super::{super::CalculatorBase, SphericalExpansionParameters};
use super::{SphericalExpansion, RadialBasis, CutoffFunction};

/// Parameters for radial spectrum calculator.
///
/// In the radial spectrum, each sample represents the spherically averaged
/// atomic density projected onto a radial basis.  In practice, this is
/// equivalent to selecting the `l=0` channel of the spherical expansion.
///
/// See [this review article](https://doi.org/10.1063/1.5090481) for more
/// information on the SOAP-like representations.
#[derive(Debug, Clone)]
#[derive(serde::Deserialize, serde::Serialize, schemars::JsonSchema)]
pub struct RadialSpectrumParameters {
    /// Spherical cutoff to use for atomic environments
    pub cutoff: f64,
    /// Number of radial basis function to use
    pub max_radial: usize,
    /// Width of the atom-centered gaussian creating the atomic density
    pub atomic_gaussian_width: f64,
    /// Should we also compute gradients of the feature?
    pub gradients: bool,
    /// radial basis to use for the radial integral
    pub radial_basis: RadialBasis,
    /// cutoff function used to smooth the behavior around the cutoff radius
    pub cutoff_function: CutoffFunction,
}

pub struct RadialSpectrum {
    parameters: RadialSpectrumParameters,
    spherical_expansion_calculator: Calculator,
    spherical_expansion: Descriptor,
}

impl RadialSpectrum {
    pub fn new(parameters: RadialSpectrumParameters) -> Result<RadialSpectrum, Error> {

        let expansion_parameters = SphericalExpansionParameters {
            cutoff: parameters.cutoff,
            max_radial: parameters.max_radial,
            max_angular: 0,
            atomic_gaussian_width: parameters.atomic_gaussian_width,
            gradients: parameters.gradients,
            radial_basis: parameters.radial_basis,
            cutoff_function: parameters.cutoff_function,
        };

        let spherical_expansion = SphericalExpansion::new(expansion_parameters)?;

        return Ok(RadialSpectrum{
            parameters: parameters,
            spherical_expansion_calculator: Calculator::from(
                Box::new(spherical_expansion) as Box<dyn CalculatorBase>),
            spherical_expansion: Descriptor::new(),
        })
    }
}

impl CalculatorBase for RadialSpectrum {
    /// Get the name of this Calculator
    fn name(&self) -> String {todo!()}

    /// Get the parameters used to create this Calculator as a JSON string
    fn get_parameters(&self) -> String {todo!()}

    /// Get the names of features for this Calculator
    fn features_names(&self) -> Vec<&str> {todo!()}
    /// Get the default set of features for this Calculator
    fn features(&self) -> Indexes {todo!()}

    /// Get the default sample builder for this Calculator
    fn samples_builder(&self) -> Box<dyn SamplesBuilder> {todo!()}
    /// Does this calculator compute gradients?
    fn compute_gradients(&self) -> bool {todo!()}

    /// Check that the given indexes are valid feature indexes for this
    /// Calculator. This is used by to ensure only valid features are requested
    fn check_features(&self, indexes: &Indexes) -> Result<(), Error> {todo!()}

    /// Check that the given indexes are valid samples indexes for this
    /// Calculator. This is used by to ensure only valid samples are requested
    ///
    /// The default implementation recompute the full set of samples using
    /// `Self::samples()`, and check that all requested samples are part of the
    /// full sample set.
    fn check_samples(&self, indexes: &Indexes, systems: &mut [Box<dyn System>]) -> Result<(), Error> {todo!()}
    fn compute(&mut self, systems: &mut [Box<dyn System>], descriptor: &mut Descriptor) -> Result<(), Error> {todo!()}
}
