#![cfg(test)]

use crate::homology::{homology_groups, HomologyGroup};
use crate::parser::parse_complex;
use crate::simplex::SimplicialComplex;
use crate::smith::smith_normal_form;
use ndarray::array;

#[test]
fn test_parser() {
    let test_cases = vec![
        ("[[1, 2], [3, 4]]", vec![vec![1, 2], vec![3, 4]]),
        (
            "[ [3, 5, 7], [9]
        ]",
            vec![vec![3, 5, 7], vec![9]],
        ),
    ];
    for (str, res) in test_cases.iter() {
        assert_eq!(&parse_complex(str).unwrap().1, res);
    }
}

#[test]
fn test_smith() {
    let test_cases = vec![
        array![[1, 2, 3,], [4, 5, 6,], [7, 8, 9],],
        array![[-1, 1, 1, 43]],
        array![
            [12, -65, 22, 0],
            [-7, 0, 43, 18],
            [68, -39, 2, 0],
            [-11, -11, -11, -11],
        ],
        array![[128], [64], [32], [16], [-8]],
    ];
    for case in test_cases.iter() {
        let (s, a, t) = smith_normal_form(case.view());
        // Check decomposition
        assert_eq!(&s.dot(&a).dot(&t), case);
        // Check non-diagonal
        a.indexed_iter().all(|((i, j), x)| i == j || *x == 0);
        // Check invariant order
        let invariants: Vec<i128> = a.diag().iter().map(|x| *x).filter(|x| *x != 0).collect();
        for i in 1..invariants.len() {
            assert_eq!(invariants[i] % invariants[i - 1], 0);
        }
    }
}

#[test]
fn test_homology() {
    let test_cases = vec![(
        // Klein bottle
        "[
            [1, 4, 6],
            [1, 2, 6],
            [2, 6, 7],
            [2, 3, 7],
            [1, 3, 7],
            [1, 4, 7],
            [4, 5, 9],
            [4, 6, 9],
            [6, 8, 9],
            [6, 7, 8],
            [5, 7, 8],
            [4, 5, 7],
            [1, 5, 9],
            [1, 3, 9],
            [2, 3, 9],
            [2, 8, 9],
            [1, 2, 8],
            [1, 5, 8] 
            ]",
        vec![
            HomologyGroup {
                free: 1,
                torsion: vec![],
            },
            HomologyGroup {
                free: 1,
                torsion: vec![2],
            },
        ],
    )];

    for (str, hom) in test_cases.iter() {
        let simp = SimplicialComplex::new(parse_complex(str).unwrap().1);
        assert_eq!(&homology_groups(simp), hom);
    }
}
