use std::collections::BTreeSet;

use crate::{Error, System};
use super::super::{SamplesBuilder, Indexes, IndexesBuilder, IndexValue};

/// `StructureSpeciesSamples` is used to represents samples corresponding to
/// full structures, where each chemical species in the structure is represented
/// separately.
///
/// The base set of indexes contains `structure` and `species` the  gradient
/// indexes also contains the `atom` inside the structure with respect to which
/// the gradient is taken and the `spatial` (i.e. x/y/z) index.
pub struct StructureSpeciesSamples;

impl SamplesBuilder for StructureSpeciesSamples {
    fn names(&self) -> Vec<&str> {
        vec!["structure", "species"]
    }

    fn gradients_names(&self) -> Option<Vec<&str>> {
        let mut names = self.names();
        names.extend_from_slice(&["atom", "spatial"]);
        return Some(names);
    }

    fn samples(&self, systems: &mut [Box<dyn System>]) -> Result<Indexes, Error> {
        let mut indexes = IndexesBuilder::new(self.names());
        for (i_system, system) in systems.iter().enumerate() {
            for &species in system.species()?.iter().collect::<BTreeSet<_>>() {
                indexes.add(&[
                    IndexValue::from(i_system), IndexValue::from(species)
                ]);
            }
        }
        return Ok(indexes.finish());
    }

    fn gradients_for(&self, systems: &mut [Box<dyn System>], samples: &Indexes) -> Result<Option<Indexes>, Error> {
        assert_eq!(samples.names(), self.names());

        let mut gradients = IndexesBuilder::new(self.gradients_names().expect("gradients names"));
        for value in samples.iter() {
            let i_system = value[0];
            let alpha = value[1];

            let system = &systems[i_system.usize()];
            let species = system.species()?;
            for (i_atom, &species) in species.iter().enumerate() {
                // only atoms with the same species participate to the gradient
                if species == alpha.i32() {
                    gradients.add(&[i_system, alpha, IndexValue::from(i_atom), IndexValue::from(0)]);
                    gradients.add(&[i_system, alpha, IndexValue::from(i_atom), IndexValue::from(1)]);
                    gradients.add(&[i_system, alpha, IndexValue::from(i_atom), IndexValue::from(2)]);
                }
            }
        }

        return Ok(Some(gradients.finish()));
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::systems::test_utils::test_systems;

    // small helper function to create IndexValue
    fn v(i: i32) -> IndexValue { IndexValue::from(i) }

    #[test]
    fn samples() {
        let mut systems = test_systems(&["methane", "methane", "water"]);
        let builder = StructureSpeciesSamples;
        assert_eq!(builder.names(), &["structure", "species"]);

        let samples = builder.samples(&mut systems).unwrap();
        assert_eq!(samples.names(), builder.names());
        assert_eq!(samples.count(), 6);
        assert_eq!(samples.iter().collect::<Vec<_>>(), vec![
            &[v(0), v(1)], &[v(0), v(6)],
            &[v(1), v(1)], &[v(1), v(6)],
            &[v(2), v(1)], &[v(2), v(123456)],
        ]);
    }

    #[test]
    fn gradients() {
        let mut systems = test_systems(&["CH", "water"]);
        let builder = StructureSpeciesSamples;
        assert_eq!(builder.gradients_names().unwrap(), &["structure", "species", "atom", "spatial"]);

        let (_, gradients) = StructureSpeciesSamples.with_gradients(&mut systems).unwrap();
        let gradients = gradients.unwrap();
        assert_eq!(gradients.count(), 15);
        assert_eq!(gradients.names(), builder.gradients_names().unwrap());

        assert_eq!(gradients.iter().collect::<Vec<_>>(), vec![
            // H channel in CH
            &[v(0), v(1), v(0), v(0)], &[v(0), v(1), v(0), v(1)], &[v(0), v(1), v(0), v(2)],
            // C channel in CH
            &[v(0), v(6), v(1), v(0)], &[v(0), v(6), v(1), v(1)], &[v(0), v(6), v(1), v(2)],
            // H channel in water
            &[v(1), v(1), v(1), v(0)], &[v(1), v(1), v(1), v(1)], &[v(1), v(1), v(1), v(2)],
            &[v(1), v(1), v(2), v(0)], &[v(1), v(1), v(2), v(1)], &[v(1), v(1), v(2), v(2)],
            // O channel in water
            &[v(1), v(123456), v(0), v(0)], &[v(1), v(123456), v(0), v(1)], &[v(1), v(123456), v(0), v(2)],
        ]);
    }

    #[test]
    fn partial_gradients() {
        let mut samples = IndexesBuilder::new(vec!["structure", "species"]);
        samples.add(&[v(2), v(1)]);
        samples.add(&[v(0), v(6)]);
        let samples = samples.finish();

        let mut systems = test_systems(&["CH", "water", "CH"]);
        let gradients = StructureSpeciesSamples.gradients_for(&mut systems, &samples).unwrap();
        let gradients = gradients.unwrap();
        assert_eq!(gradients.names(), StructureSpeciesSamples.gradients_names().unwrap());

        assert_eq!(gradients.iter().collect::<Vec<_>>(), vec![
            // H channel in CH #2
            &[v(2), v(1), v(0), v(0)],
            &[v(2), v(1), v(0), v(1)],
            &[v(2), v(1), v(0), v(2)],
            // C channel in CH #1
            &[v(0), v(6), v(1), v(0)],
            &[v(0), v(6), v(1), v(1)],
            &[v(0), v(6), v(1), v(2)],
        ]);
    }
}
