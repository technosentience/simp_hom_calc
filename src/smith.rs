use ndarray::{s, Array2, ArrayView1, ArrayView2, ArrayViewMut2, Zip};
use std::cmp::min;

fn extended_euclid(a: i128, b: i128) -> (i128, i128, i128) {
    if a == 0 {
        return (b, 0, 1);
    } else {
        let (gcd, x, y) = extended_euclid(b.rem_euclid(a), a);
        return (gcd, y - b.div_euclid(a) * x, x);
    }
}

// Left-multuply by matrix with 2 nonzero rows/cols
fn left_mul2<'a>(
    mut m: ArrayViewMut2<'a, i128>,
    i: usize,
    j: usize,
    a: i128,
    b: i128,
    c: i128,
    d: i128,
) {
    let (row_i, row_j) = m.multi_slice_mut((s![i, ..], s![j, ..]));
    Zip::from(row_i).and(row_j).for_each(|x, y| {
        let (xval, yval) = (*x, *y);
        *x = a * xval + b * yval;
        *y = c * xval + d * yval;
    });
}

// Right-multuply by matrix with 2 nonzero rows/cols
fn right_mul2<'a>(
    mut m: ArrayViewMut2<'a, i128>,
    i: usize,
    j: usize,
    a: i128,
    b: i128,
    c: i128,
    d: i128,
) {
    let (col_i, col_j) = m.multi_slice_mut((s![.., i], s![.., j]));
    Zip::from(col_i).and(col_j).for_each(|x, y| {
        let (xval, yval) = (*x, *y);
        *x = a * xval + c * yval;
        *y = b * xval + d * yval;
    });
}

fn find_nonzero_entry<'a>(v: ArrayView1<'a, i128>) -> Option<usize> {
    v.iter().position(|x| *x != 0)
}

fn find_nonzero_column<'a>(m: ArrayView2<'a, i128>) -> Option<usize> {
    m.columns()
        .into_iter()
        .position(|c| find_nonzero_entry(c) != None)
}

fn make_pivot<'a>(
    mut s: ArrayViewMut2<'a, i128>,
    mut a: ArrayViewMut2<'a, i128>,
    mut t: ArrayViewMut2<'a, i128>,
    i: usize,
    start_j: usize,
) -> Option<usize> {
    let j = find_nonzero_column(a.slice(s![.., start_j..]))? + start_j;

    let nonzero_i = find_nonzero_entry(a.column(j))?;
    if i != nonzero_i {
        // Swap rows
        left_mul2(a.view_mut(), i, nonzero_i, 0, 1, 1, 0);
        right_mul2(s.view_mut(), i, nonzero_i, 0, 1, 1, 0);
    }
    if i != j {
        // Swap columns
        right_mul2(a.view_mut(), i, j, 0, 1, 1, 0);
        left_mul2(t.view_mut(), i, j, 0, 1, 1, 0);
    }
    Some(j + 1)
}

fn zero_row_entry<'a>(
    mut a: ArrayViewMut2<'a, i128>,
    mut s: ArrayViewMut2<'a, i128>,
    i0: usize,
    i1: usize,
) {
    let (gcd, x, y) = extended_euclid(a[(i0, i0)], a[(i1, i0)]);
    let z = a[(i1, i0)] / gcd;
    let w = a[(i0, i0)] / gcd;
    // a[i0, i0] = gcd, a[i1, i0] = 0
    left_mul2(a.view_mut(), i0, i1, x, y, -z, w);
    right_mul2(s.view_mut(), i0, i1, w, -y, z, x);
}

fn zero_column_entry<'a>(
    mut a: ArrayViewMut2<'a, i128>,
    mut t: ArrayViewMut2<'a, i128>,
    i0: usize,
    i1: usize,
) {
    let (gcd, x, y) = extended_euclid(a[(i0, i0)], a[(i0, i1)]);
    let z = a[(i0, i1)] / gcd;
    let w = a[(i0, i0)] / gcd;
    // a[i0, i0] = gcd, a[i0, i1] = 0
    right_mul2(a.view_mut(), i0, i1, x, -z, y, w);
    left_mul2(t.view_mut(), i0, i1, w, z, -y, x);
}

fn reduce_row<'a>(
    mut s: ArrayViewMut2<'a, i128>,
    mut a: ArrayViewMut2<'a, i128>,
    mut t: ArrayViewMut2<'a, i128>,
    i0: usize,
    start_j: usize,
) -> Option<usize> {
    let next_j = make_pivot(s.view_mut(), a.view_mut(), t.view_mut(), i0, start_j)?;
    loop {
        if let Some(i1) = find_nonzero_entry(a.slice(s![i0 + 1.., i0])) {
            zero_row_entry(a.view_mut(), s.view_mut(), i0, i1 + i0 + 1);
        } else if let Some(i1) = find_nonzero_entry(a.slice(s![i0, i0 + 1..])) {
            zero_column_entry(a.view_mut(), t.view_mut(), i0, i1 + i0 + 1);
        } else {
            break;
        }
    }
    Some(next_j)
}

fn reduce_diagonal<'a>(
    mut s: ArrayViewMut2<'a, i128>,
    mut a: ArrayViewMut2<'a, i128>,
    mut t: ArrayViewMut2<'a, i128>,
    i1: usize,
) {
    for i0 in 0..i1 {
        let (gcd, x, y) = extended_euclid(a[(i0, i0)], a[(i1, i1)]);
        if gcd == 0 {
            continue;
        }
        let z = a[(i1, i1)] / gcd;
        let w = a[(i0, i0)] / gcd;
        // a[i0, i0] = gcd, a[i1, i1] = lcm
        left_mul2(a.view_mut(), i0, i1, 1, y, z, y * z - 1);
        right_mul2(s.view_mut(), i0, i1, 1 - y * z, y, z, -1);
        right_mul2(a.view_mut(), i0, i1, x, 1 - x * w, 1, -w);
        left_mul2(t.view_mut(), i0, i1, w, 1 - x * w, 1, -x);
    }
}

fn smith_normal_form_in_place<'a, 'b>(
    mut s: ArrayViewMut2<'a, i128>,
    mut a: ArrayViewMut2<'a, i128>,
    mut t: ArrayViewMut2<'a, i128>,
) {
    let (rows, cols) = a.dim();
    let mut j = 0;
    for i in 0..rows {
        let next_j = reduce_row(s.view_mut(), a.view_mut(), t.view_mut(), i, j);
        if let Some(next_j) = next_j {
            j = next_j;
        }
    }
    for i in 1..min(rows, cols) {
        reduce_diagonal(s.view_mut(), a.view_mut(), t.view_mut(), i);
    }

    for i in 0..min(rows, cols) {
        if a[(i, i)] < 0 {
            *&mut (a.row_mut(i)) *= -1;
            *&mut (s.column_mut(i)) *= -1;
        }
    }
}

/// Compute the Smith normal form M = S A T, where all matrices are integer,
/// S and T are invertible and A has only diagonal entries with a_i | a_i+1.
pub fn smith_normal_form(m: ArrayView2<i128>) -> (Array2<i128>, Array2<i128>, Array2<i128>) {
    let mut a = m.into_owned();
    let mut s = Array2::eye(m.nrows());
    let mut t = Array2::eye(m.ncols());
    smith_normal_form_in_place(s.view_mut(), a.view_mut(), t.view_mut());
    return (s, a, t);
}
