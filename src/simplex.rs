use bitvec::prelude::*;
use ndarray::Array2;
use std::collections::{HashMap, HashSet, VecDeque};

pub struct SimplicialComplex {
    simplices: HashSet<BitVec>,
}

impl SimplicialComplex {
    pub fn dimension(&self) -> usize {
        return self.simplices.iter().map(|s| s.count_ones()).max().unwrap() - 1;
    }

    fn simplex_of_vec(points: &Vec<usize>, bit_size: usize) -> BitVec {
        let mut simp = BitVec::new();
        simp.resize(bit_size, false);
        for i in points {
            simp.set(*i, true);
        }
        return simp;
    }

    fn iter_over_subsets<F>(simp: &BitVec, mut act: F)
    where
        F: FnMut(&BitVec, bool),
    {
        let mut subset = BitVec::new();
        let mut sign = true;
        subset.resize(simp.len(), false);
        subset.copy_from_bitslice(&simp);
        for i in 0..subset.len() {
            if subset[i] {
                subset.set(i, false);
                act(&subset, sign);
                subset.set(i, true);
                sign = !sign;
            }
        }
    }

    pub fn boundary_map(&self, dim: usize) -> Array2<i128> {
        let target: HashMap<&BitVec, usize> = self
            .simplices
            .iter()
            .filter(|s| s.count_ones() == dim)
            .enumerate()
            .map(|(i, s)| (s, i))
            .collect();

        let source: Vec<&BitVec> = self
            .simplices
            .iter()
            .filter(|s| s.count_ones() == dim + 1)
            .collect();

        let mut map = Array2::<i128>::zeros((target.len(), source.len()));
        for (j, s) in source.iter().enumerate() {
            Self::iter_over_subsets(s, |subset, sign| {
                let i = target[subset];
                map[(i, j)] = if sign { 1 } else { -1 };
            });
        }

        return map;
    }

    pub fn new(simp_vec: Vec<Vec<usize>>) -> Self {
        let points = simp_vec
            .iter()
            .map(|s| s.iter().max().unwrap())
            .max()
            .unwrap()
            .to_owned()
            + 1;
        let mut simplices = HashSet::new();
        let mut queue: VecDeque<BitVec> = simp_vec
            .iter()
            .map(|s| Self::simplex_of_vec(s, points))
            .collect();
        while !queue.is_empty() {
            let simp = queue.pop_front().unwrap();
            if simplices.contains(&simp) {
                continue;
            }
            Self::iter_over_subsets(&simp, |sub_simp, _| {
                queue.push_back(sub_simp.to_owned());
            });
            simplices.insert(simp);
        }

        return Self { simplices };
    }
}
