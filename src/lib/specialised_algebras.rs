// mod libs {
//     include!("../lib/algebras.rs");
// }

// use crate::libs::*;

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

// N-dimensional matrix addition group
#[derive(Debug, Clone, PartialEq)]
pub struct NMatrixAddition {
    pub elements: Vec<f64>,
}

impl NMatrixAddition {
    pub fn zero() -> Self {
        NMatrixAddition { elements: vec![0.0; 2] }
    }
}

impl Algebra<NMatrixAddition> for NMatrixAddition {
    fn operate_binary(&self, a: NMatrixAddition, b: NMatrixAddition) -> NMatrixAddition {
        NMatrixAddition {
            elements: a.elements.iter().zip(b.elements.iter()).map(|(x, y)| x + y).collect(),
        }
    }

    fn operate(&self, elements: &[NMatrixAddition]) -> NMatrixAddition {
        if elements.is_empty() {
            NMatrixAddition::zero()
        } else {
            elements.iter().fold(NMatrixAddition::zero(), |acc, x| self.operate_binary(acc, (*x).clone()))
        }
    }
}
impl Group<NMatrixAddition> for NMatrixAddition {
    fn identity(&self) -> NMatrixAddition {
        NMatrixAddition::zero()
    }
    fn inverse(&self, a: NMatrixAddition) -> NMatrixAddition {
        NMatrixAddition {
            elements: a.elements.iter().map(|&x| -x).collect(),
        }
    }
}

// N-dimensional matrix multiplication group
#[derive(Debug, Clone, PartialEq)]
pub struct NMatrixMultiplication {
    pub elements: Vec<f64>,
}

impl Algebra<NMatrixMultiplication> for NMatrixMultiplication {
    fn operate_binary(&self, a: NMatrixMultiplication, b: NMatrixMultiplication) -> NMatrixMultiplication {
        NMatrixMultiplication {
            elements: a.elements.iter().zip(b.elements.iter()).map(|(x, y)| x * y).collect(),
        }
    }

    fn operate(&self, elements: &[NMatrixMultiplication]) -> NMatrixMultiplication {
        if elements.is_empty() {
            NMatrixMultiplication::zero()
        } else {
            elements.iter().fold(NMatrixMultiplication::zero(), |acc, x| self.operate_binary(acc, (*x).clone()))
        }
    }
}

impl NMatrixMultiplication {
    pub fn zero() -> Self {
        NMatrixMultiplication { elements: vec![1.0; 2] } // Identity for multiplication
    }
}

impl Group<NMatrixMultiplication> for NMatrixMultiplication {
    fn identity(&self) -> NMatrixMultiplication {
        NMatrixMultiplication::zero()
    }
    fn inverse(&self, a: NMatrixMultiplication) -> NMatrixMultiplication {
        NMatrixMultiplication {
            elements: a.elements.iter().map(|&x| if x != 0.0 { 1.0 / x } else { 0.0 }).collect(),
        }
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

#[derive(Debug, Clone, Copy, PartialEq)]
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

impl Algebra<Octonion> for OctonionMultiplication {
    fn operate_binary(&self, a: Octonion, b: Octonion) -> Octonion {
        Octonion {
            re: a.re * b.re - a.i * b.i - a.j * b.j - a.k * b.k - a.l * b.l - a.m * b.m - a.n * b.n - a.o * b.o,
            i: a.re * b.i + a.i * b.re,
            j: a.re * b.j + a.j * b.re,
            k: a.re * b.k + a.k * b.re,
            l: a.re * b.l + a.l * b.re,
            m: a.re * b.m + a.m * b.re,
            n: a.re * b.n + a.n * b.re,
            o: a.re * b.o + a.o * b.re,
        }
    }

    fn operate(&self, elements: &[Octonion]) -> Octonion {
        elements.iter().fold(self.identity(), |acc, &x| self.operate_binary(acc, x))
    }
}

impl Group<Octonion> for OctonionMultiplication {
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

    fn inverse(&self, a: Octonion) -> Octonion {
        let norm_squared = a.re * a.re + a.i * a.i + a.j * a.j + a.k * a.k + a.l * a.l + a.m * a.m + a.n * a.n + a.o * a.o;
        if norm_squared == 0.0 {
            panic!("Cannot invert a zero octonion");
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
}