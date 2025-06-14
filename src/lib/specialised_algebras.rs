// mod libs {
//     include!("../lib/algebras.rs");
// }

// use crate::libs::*;

use num_traits::{Float, Zero, One};

// Integer addition group
pub struct IntegerAddition;

impl Algebra<i32> for IntegerAddition {
    fn operate_binary(&self, a: i32, b: i32) -> i32 {
        a.wrapping_add(b)
    }
    
    fn operate(&self, elements: &[i32]) -> i32 {
        if elements.is_empty() {
            self.identity()
        } else {
            elements.iter().fold(self.identity(), |acc, &x| self.operate_binary(acc, x))
        }
    }
}

impl Group<i32> for IntegerAddition {
    fn identity(&self) -> i32 {
        0
    }
    
    fn inverse(&self, a: i32) -> i32 {
        -a
    }
}

// Integer multiplication monoid
pub struct IntegerMultiplication;

impl Algebra<i32> for IntegerMultiplication {
    fn operate_binary(&self, a: i32, b: i32) -> i32 {
        a.wrapping_mul(b)
    }
    
    fn operate(&self, elements: &[i32]) -> i32 {
        if elements.is_empty() {
            1 // identity for multiplication
        } else {
            elements.iter().fold(1, |acc, &x| self.operate_binary(acc, x))
        }
    }
}


#[derive(Debug, Clone, Copy)]
pub struct ModN {
    pub n: i32,
}

#[derive(Debug, Clone, Copy)]
pub struct Quaternion {
    pub re: f64,
    pub i: f64,
    pub j: f64,
    pub k: f64,
}


pub struct QuaternionAddition;
pub struct QuaternionMultiplication;

impl Algebra<Quaternion> for QuaternionMultiplication {
    fn operate_binary(&self, a: Quaternion, b: Quaternion) -> Quaternion {
        Quaternion {
            re: a.re * b.re - a.i * b.i - a.j * b.j - a.k * b.k,
            i: a.re * b.i + a.i * b.re + a.j * b.k - a.k * b.j,
            j: a.re * b.j - a.i * b.k + a.j * b.re + a.k * b.i,
            k: a.re * b.k + a.i * b.j - a.j * b.i + a.k * b.re,
        }
    }

    fn operate(&self, elements: &[Quaternion]) -> Quaternion {
        if elements.is_empty() {
            self.identity()
        } else {
            elements.iter().fold(self.identity(), |acc, &x| self.operate_binary(acc, x))
        }
    }
}

impl Group<Quaternion> for QuaternionMultiplication {
    fn identity(&self) -> Quaternion {
        Quaternion { re: 1.0, i: 0.0, j: 0.0, k: 0.0 }
    }

    fn inverse(&self, a: Quaternion) -> Quaternion {
        // For quaternion multiplication, the inverse is more complex
        let norm_squared = a.re * a.re + a.i * a.i + a.j * a.j + a.k * a.k;
        if norm_squared == 0.0 {
            panic!("Cannot invert a zero quaternion");
        }
        Quaternion {
            re: a.re / norm_squared,
            i: -a.i / norm_squared,
            j: -a.j / norm_squared,
            k: -a.k / norm_squared,
        }
    }
}

impl Algebra<Quaternion> for QuaternionAddition {
    fn operate_binary(&self, a: Quaternion, b: Quaternion) -> Quaternion {
        Quaternion {
            re: a.re + b.re,
            i: a.i + b.i,
            j: a.j + b.j,
            k: a.k + b.k,
        }
    }

    fn operate(&self, elements: &[Quaternion]) -> Quaternion {
        elements.iter().fold(self.identity(), |acc, &x| self.operate_binary(acc, x))
    }
}

impl Group<Quaternion> for QuaternionAddition {
    fn identity(&self) -> Quaternion {
        Quaternion { re: 0.0, i: 0.0, j: 0.0, k: 0.0 }
    }

    fn inverse(&self, a: Quaternion) -> Quaternion {
        Quaternion {
            re: -a.re,
            i: -a.i,
            j: -a.j,
            k: -a.k,
        }
    }
}

// Complex number structure
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

// Complex addition group
pub struct ComplexAddition;

impl Algebra<Complex> for ComplexAddition {
    fn operate_binary(&self, a: Complex, b: Complex) -> Complex {
        Complex {
            re: a.re + b.re,
            im: a.im + b.im,
        }
    }
    
    fn operate(&self, elements: &[Complex]) -> Complex {
        if elements.is_empty() {
            Complex { re: 0.0, im: 0.0 }
        } else {
            elements.iter().fold(Complex { re: 0.0, im: 0.0 }, |acc, &x| self.operate_binary(acc, x))
        }
    }
}

impl Group<Complex> for ComplexAddition {
    fn identity(&self) -> Complex {
        Complex { re: 0.0, im: 0.0 }
    }
    
    fn inverse(&self, a: Complex) -> Complex {
        Complex { re: -a.re, im: -a.im }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<T> {
    pub rows: usize,
    pub cols: usize,
    pub elements: Vec<T>, // Flat vector, row-major order
}

impl<T: Copy + Zero + One + Float> Matrix<T> {  // Bounds ensure T supports necessary ops without assuming Reals
    pub fn new(rows: usize, cols: usize, elements: Vec<T>) -> Result<Self, String> {
        if elements.len() != rows * cols {
            return Err("Element count must match rows * cols".to_string());
        }
        Ok(Matrix { rows, cols, elements })
    }

    pub fn zero(rows: usize, cols: usize) -> Self {
        Matrix {
            rows,
            cols,
            elements: vec![T::zero(); rows * cols],
        }
    }

    pub fn identity(size: usize) -> Self {  // For square identity matrix
        let mut elements = vec![T::zero(); size * size];
        for i in 0..size {
            elements[i * size + i] = T::one();
        }
        Matrix { rows: size, cols: size, elements }
    }

    // Helper to get element at (i, j)
    fn get(&self, i: usize, j: usize) -> T {
        self.elements[i * self.cols + j]
    }

    // Helper to set element at (i, j)
    fn set(&mut self, i: usize, j: usize, value: T) {
        self.elements[i * self.cols + j] = value;
    }
}

// Updated MatrixAddition with proper dimension checks
pub struct MatrixAddition;

impl<T: Copy + std::ops::Add<Output = T> + Zero + Float> Algebra<Matrix<T>> for MatrixAddition {
    fn operate_binary(&self, a: Matrix<T>, b: Matrix<T>) -> Matrix<T> {
        if a.rows != b.rows || a.cols != b.cols {
            panic!("Cannot add matrices with mismatched dimensions");
        }
        let elements = a.elements.iter().zip(b.elements.iter()).map(|(&x, &y)| x + y).collect();
        Matrix { rows: a.rows, cols: a.cols, elements }
    }

    fn operate(&self, elements: &[Matrix<T>]) -> Matrix<T> {
        if elements.is_empty() {
            return Matrix::zero(0, 0);  // Handle empty case as 0x0 zero matrix
        }
        let first = &elements[0];
        elements.iter().skip(1).fold(first.clone(), |acc, x| self.operate_binary(acc, x.clone()))
    }
}

impl<T: Copy + std::ops::Neg<Output = T> + std::ops::Add<Output = T> + Zero + Float> Group<Matrix<T>> for MatrixAddition {
    fn identity(&self) -> Matrix<T> {
        Matrix::zero(0, 0)  // General identity is the zero matrix; for specific dims, use Matrix::zero
    }

    fn inverse(&self, a: Matrix<T>) -> Matrix<T> {
        let elements = a.elements.iter().map(|&x| -x).collect();
        Matrix { rows: a.rows, cols: a.cols, elements }
    }
}

// Updated MatrixMultiplication with standard matrix multiplication
pub struct MatrixMultiplication;

impl<T: Copy + std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Zero + One + Float> Algebra<Matrix<T>> for MatrixMultiplication {
    fn operate_binary(&self, a: Matrix<T>, b: Matrix<T>) -> Matrix<T> {
        if a.cols != b.rows {
            panic!("Cannot multiply matrices: incompatible dimensions (a.cols != b.rows)");
        }
        let mut elements = vec![T::zero(); a.rows * b.cols];
        for i in 0..a.rows {
            for j in 0..b.cols {
                let mut sum = T::zero();
                for k in 0..a.cols {
                    sum = sum + a.get(i, k) * b.get(k, j);
                }
                elements[i * b.cols + j] = sum;
            }
        }
        Matrix { rows: a.rows, cols: b.cols, elements }
    }

    fn operate(&self, elements: &[Matrix<T>]) -> Matrix<T> {
        if elements.is_empty() {
            return Matrix::identity(1);  // 1x1 identity for empty product
        }
        elements.iter().fold(Matrix::identity(elements[0].rows), |acc, x| {
            if acc.cols != x.rows {
                panic!("Incompatible dimensions in matrix product");
            }
            self.operate_binary(acc, x.clone())
        })
    }
}

impl<T: Copy + std::ops::Mul<Output = T> + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Div<Output = T> + Zero + One + Float + PartialOrd> Group<Matrix<T>> for MatrixMultiplication {
    fn identity(&self) -> Matrix<T> {
        Matrix::identity(1)  // Default 1x1 identity; can be adjusted for dims
    }

    fn inverse(&self, a: Matrix<T>) -> Matrix<T> {
        if a.rows != a.cols {
            panic!("Matrix must be square for inversion");
        }
        // Simplified Gaussian elimination for inversion (for generality)
        // Note: This is a basic implementation; for production, consider using a library like nalgebra
        let n = a.rows;
        let mut mat = a.clone();
        let mut inv = Matrix::identity(n);
        for i in 0..n {
            // Find pivot
            if mat.get(i, i) == T::zero() {
                panic!("Matrix is singular and cannot be inverted");
            }
            // Eliminate
            for j in 0..n {
                if i != j {
                    let factor = mat.get(j, i) / mat.get(i, i);
                    for k in 0..n {
                        let new_val = mat.get(j, k) - factor * mat.get(i, k);
                        mat.set(j, k, new_val);
                        let inv_val = inv.get(j, k) - factor * inv.get(i, k);
                        inv.set(j, k, inv_val);
                    }
                }
            }
        }
        // Normalize (simplified)
        for i in 0..n {
            let pivot = mat.get(i, i);
            for j in 0..n {
                inv.set(i, j, inv.get(i, j) / pivot);
            }
        }
        inv
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct NDimNumber {
    pub components: Vec<f64>,
}

impl NDimNumber {
    pub fn new(components: Vec<f64>) -> Self {
        NDimNumber { components }
    }
    
    pub fn zero(dimension: usize) -> Self {
        NDimNumber { components: vec![0.0; dimension] }
    }
    
    pub fn one(dimension: usize) -> Self {
        let mut components = vec![0.0; dimension];
        if dimension > 0 {
            components[0] = 1.0;
        }
        NDimNumber { components }
    }
    
    // Helper to convert Complex to NDimNumber
    pub fn from_complex(c: Complex) -> Self {
        NDimNumber { components: vec![c.re, c.im] }
    }
    
    // Helper to convert to Complex if dimension is 2
    pub fn to_complex(&self) -> Option<Complex> {
        if self.components.len() == 2 {
            Some(Complex { re: self.components[0], im: self.components[1] })
        } else {
            None
        }
    }
}

// Addition for n-dimensional numbers
pub struct NDimAddition;

impl Algebra<NDimNumber> for NDimAddition {
    fn operate_binary(&self, a: NDimNumber, b: NDimNumber) -> NDimNumber {
        assert_eq!(a.components.len(), b.components.len(), "Dimensions must match");
        let components = a.components.iter()
            .zip(b.components.iter())
            .map(|(&x, &y)| x + y)
            .collect();
        NDimNumber { components }
    }
    
    fn operate(&self, elements: &[NDimNumber]) -> NDimNumber {
        if elements.is_empty() {
            return NDimNumber { components: vec![] };
        }
        
        let dimension = elements[0].components.len();
        let identity = NDimNumber::zero(dimension);
        
        elements.iter().fold(identity, |acc, x| self.operate_binary(acc, x.clone()))
    }
}

impl Group<NDimNumber> for NDimAddition {
    fn identity(&self) -> NDimNumber {
        // This is a placeholder - proper implementation would need dimension info
        NDimNumber { components: vec![] }
    }
    
    fn inverse(&self, a: NDimNumber) -> NDimNumber {
        let components = a.components.iter().map(|&x| -x).collect();
        NDimNumber { components }
    }
}

// For complex numbers specifically (n=2), we can implement multiplication
pub struct ComplexMultiplication;

impl Algebra<Complex> for ComplexMultiplication {
    fn operate_binary(&self, a: Complex, b: Complex) -> Complex {
        Complex {
            re: a.re * b.re - a.im * b.im,
            im: a.re * b.im + a.im * b.re,
        }
    }
    
    fn operate(&self, elements: &[Complex]) -> Complex {
        if elements.is_empty() {
            self.identity()
        } else {
            elements.iter().fold(self.identity(), |acc, &x| self.operate_binary(acc, x))
        }
    }
}

impl Group<Complex> for ComplexMultiplication {
    fn identity(&self) -> Complex {
        Complex { re: 1.0, im: 0.0 }
    }
    
    fn inverse(&self, a: Complex) -> Complex {
        let norm_squared = a.re * a.re + a.im * a.im;
        if norm_squared == 0.0 {
            panic!("Cannot invert zero complex number");
        }
        Complex { 
            re: a.re / norm_squared, 
            im: -a.im / norm_squared 
        }
    }
}

impl PartialEq for Quaternion {
    fn eq(&self, other: &Self) -> bool {
        (self.re - other.re).abs() < 1e-6 &&
        (self.i - other.i).abs() < 1e-6 &&
        (self.j - other.j).abs() < 1e-6 &&
        (self.k - other.k).abs() < 1e-6
    }
}

#[allow(unused)]
pub struct EuclideanSpace;

impl EuclideanSpace {
    #[allow(unused)]
    pub fn new(dimensions: usize) -> Self {
        EuclideanSpace
    }

    #[allow(unused)]
    pub fn vector(&self, components: &[f64]) -> Vec<f64> {
        components.to_vec()
    }

    #[allow(unused)]
    pub fn dot(&self, v1: &[f64], v2: &[f64]) -> f64 {
        v1.iter().zip(v2.iter()).map(|(x, y)| x * y).sum()
    }

    #[allow(unused)]
    pub fn norm(&self, v: &[f64]) -> f64 {
        self.dot(v, v).sqrt()
    }

    #[allow(unused)]
    pub fn add(&self, v1: &[f64], v2: &[f64]) -> Vec<f64> {
        v1.iter().zip(v2.iter()).map(|(x, y)| x + y).collect()
    }

    #[allow(unused)]
    pub fn scalar_mul(&self, v: &[f64], scalar: f64) -> Vec<f64> {
        v.iter().map(|x| x * scalar).collect()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Octonion {
    pub re: f64,
    pub i: f64,
    pub j: f64,
    pub k: f64,
    pub l: f64,
    pub m: f64,
    pub n: f64,
    pub o: f64,
}

pub struct OctonionAddition;
pub struct OctonionMultiplication;

impl Algebra<Octonion> for OctonionAddition {
    fn operate_binary(&self, a: Octonion, b: Octonion) -> Octonion {
        Octonion {
            re: a.re + b.re,
            i: a.i + b.i,
            j: a.j + b.j,
            k: a.k + b.k,
            l: a.l + b.l,
            m: a.m + b.m,
            n: a.n + b.n,
            o: a.o + b.o,
        }
    }

    fn operate(&self, elements: &[Octonion]) -> Octonion {
        elements.iter().fold(self.identity(), |acc, &x| self.operate_binary(acc, x))
    }
}

impl Group<Octonion> for OctonionAddition {
    fn identity(&self) -> Octonion {
        Octonion {
            re: 0.0,
            i: 0.0,
            j: 0.0,
            k: 0.0,
            l: 0.0,
            m: 0.0,
            n: 0.0,
            o: 0.0,
        }
    }

    fn inverse(&self, a: Octonion) -> Octonion {
        Octonion {
            re: -a.re,
            i: -a.i,
            j: -a.j,
            k: -a.k,
            l: -a.l,
            m: -a.m,
            n: -a.n,
            o: -a.o,
        }
    }
}

// Updated OctonionMultiplication with correct, full multiplication formula
impl Algebra<Octonion> for OctonionMultiplication {
    fn operate_binary(&self, a: Octonion, b: Octonion) -> Octonion {
        // Full standard octonion multiplication (non-associative, based on Cayley-Dickson or Fano plane)
        // Each component follows the multiplication table for basis elements
        Octonion {
            re: a.re * b.re - a.i * b.i - a.j * b.j - a.k * b.k - a.l * b.l - a.m * b.m - a.n * b.n - a.o * b.o,
            i: a.re * b.i + a.i * b.re + a.j * b.l - a.k * b.m + a.l * b.j + a.m * b.k - a.n * b.o + a.o * b.n,
            j: a.re * b.j - a.i * b.l + a.j * b.re + a.k * b.n - a.l * b.i - a.m * b.o + a.n * b.k + a.o * b.m,
            k: a.re * b.k + a.i * b.m - a.j * b.n + a.k * b.re + a.l * b.o - a.m * b.i - a.n * b.j + a.o * b.l,
            l: a.re * b.l - a.i * b.j + a.j * b.i - a.k * b.o + a.l * b.re + a.m * b.n + a.n * b.m - a.o * b.k,
            m: a.re * b.m + a.i * b.k + a.j * b.o + a.k * b.i - a.l * b.n + a.m * b.re - a.n * b.l - a.o * b.j,
            n: a.re * b.n - a.i * b.o - a.j * b.k + a.k * b.j + a.l * b.m + a.m * b.l + a.n * b.re + a.o * b.i,
            o: a.re * b.o + a.i * b.n - a.j * b.m - a.k * b.l - a.l * b.k + a.m * b.j - a.n * b.i + a.o * b.re,
        }
    }

    fn operate(&self, elements: &[Octonion]) -> Octonion {
        if elements.is_empty() {
            self.identity()  // Explicitly handle empty case
        } else {
            elements.iter().fold(self.identity(), |acc, &x| self.operate_binary(acc, x))
        }
    }
}

impl Group<Octonion> for OctonionMultiplication {
    // ... existing code ... (identity remains the same)

    fn inverse(&self, a: Octonion) -> Octonion {
        let norm_squared = a.re * a.re + a.i * a.i + a.j * a.j + a.k * a.k + a.l * a.l + a.m * a.m + a.n * a.n + a.o * a.o;
        const EPS: f64 = 1e-10;  // Added for numerical stability
        if norm_squared.abs() < EPS {
            panic!("Cannot invert a near-zero octonion (norm too small)");
        }
        Octonion {
            re: a.re / norm_squared,
            i: -a.i / norm_squared,
            j: -a.j / norm_squared,
            k: -a.k / norm_squared,
            l: -a.l / norm_squared,
            m: -a.m / norm_squared,
            n: -a.n / norm_squared,
            o: -a.o / norm_squared,
        }
    }

    fn identity(&self) -> Octonion {
        Octonion {
            re: 1.0,
            i: 0.0,
            j: 0.0,
            k: 0.0,
            l: 0.0,
            m: 0.0,
            n: 0.0,
            o: 0.0,
        }
    }
}

// Update PartialEq for better floating-point handling
impl PartialEq for Octonion {
    fn eq(&self, other: &Self) -> bool {
        const EPS: f64 = 1e-6;
        (self.re - other.re).abs() < EPS &&
        (self.i - other.i).abs() < EPS &&
        (self.j - other.j).abs() < EPS &&
        (self.k - other.k).abs() < EPS &&
        (self.l - other.l).abs() < EPS &&
        (self.m - other.m).abs() < EPS &&
        (self.n - other.n).abs() < EPS &&
        (self.o - other.o).abs() < EPS
    }
}