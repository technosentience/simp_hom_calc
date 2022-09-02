use crate::simplex::SimplicialComplex;
use crate::smith::smith_normal_form;
use ndarray::ArrayView2;
use std::fmt::{Display, Formatter, Result};

#[derive(Debug, PartialEq, Eq)]
pub struct HomologyGroup {
    pub free: usize,
    pub torsion: Vec<usize>,
}

impl HomologyGroup {
    fn of_boundary_maps<'a>(d0: ArrayView2<'a, i128>, d1: ArrayView2<'a, i128>) -> Self {
        let (_, d0, _) = smith_normal_form(d0);
        let (_, d1, _) = smith_normal_form(d1);
        let invariants: Vec<i128> = d1.diag().iter().map(|x| *x).filter(|x| *x != 0).collect();
        let rank1 = invariants.len();
        let rank0 = d0.diag().iter().filter(|x| **x != 0).count();
        let dim = d0.dim().1;
        let torsion: Vec<usize> = invariants
            .iter()
            .map(|x| *x as usize)
            .filter(|x| *x != 1)
            .collect();
        let free = dim - rank0 - rank1;
        return Self { free, torsion };
    }

}

impl Display for HomologyGroup {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let mut str = format!("Z^{}", self.free);
        for qout in self.torsion.iter() {
            str += &format!(" x Z/{}Z", qout);
        }
        return write!(f, "{}", str);
    }
}

pub fn homology_groups(simplex: SimplicialComplex) -> Vec<HomologyGroup> {
    let mut groups = Vec::new();
    let mut d0 = simplex.boundary_map(0);
    for dim in 0..simplex.dimension() {
        let d1 = simplex.boundary_map(dim + 1);
        let hom = HomologyGroup::of_boundary_maps(d0.view(), d1.view());
        d0 = d1;
        groups.push(hom);
    }
    return groups;
}
