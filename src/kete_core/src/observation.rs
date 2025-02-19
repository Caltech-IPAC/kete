// use nalgebra::Vector3;

// use crate::state::State;

// pub trait Observation {
//     fn error(&self, state: &State) -> f64;

//     fn observer(&self) -> &State;

//     fn jd(&self) -> f64;
// }

// pub struct ObservationVec {
//     observer: State,

//     vec: [f64; 3],

//     uncertainty: Option<[f64; 3]>,
// }

// impl Observation for ObservationVec {
//     fn error(&self, state: &State) -> f64 {
//         let unc: Vector3<f64> = self.uncertainty.clone().unwrap_or([0.0; 3]).into();
//         let mut err = 0.0;

//         let vec: Vector3<f64> = self.vec.clone().into();

//         let obs_pos: Vector3<f64> = self.observer.pos.clone().into();
//         let diff: Vector3<f64> = state.pos.clone() - obs_pos;
//         diff.angle(&vec)
//     }

//     fn jd(&self) -> f64 {
//         0.0
//     }

//     fn observer(&self) -> &State {
//         &self.observer
//     }
// }
