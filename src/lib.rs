// Trait for a general algebraic structure with a binary operation
pub trait Algebra<T> {
    fn operate_binary(&self, a: T, b: T) -> T;
    fn operate(&self, elements: &[T]) -> T;
}

// Trait for a group, which extends Algebra
pub trait Group<T>: Algebra<T> {
    fn identity(&self) -> T;
    fn inverse(&self, a: T) -> T;
}

// Concrete implementation: Integer addition group
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

impl ModN {
    pub fn new(n: i32) -> Self {
        assert!(n > 0, "Modulus must be positive and non-zero");
        ModN { n }
    }
}

impl Algebra<i32> for ModN {
    fn operate_binary(&self, a: i32, b: i32) -> i32 {
        ((a % self.n + b % self.n) % self.n + self.n) % self.n
    }

    fn operate(&self, elements: &[i32]) -> i32 {
        elements.iter().fold(self.identity(), |acc, &x| self.operate_binary(acc, x))
    }
}

impl Group<i32> for ModN {
    fn identity(&self) -> i32 {
        0
    }

    fn inverse(&self, a: i32) -> i32 {
        ((self.n - (a % self.n)) % self.n + self.n) % self.n
    }
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

// 2x2 Matrix structure
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix2x2(pub [[i32; 2]; 2]);

// Matrix addition group
pub struct MatrixAddition;

impl Algebra<Matrix2x2> for MatrixAddition {
    fn operate_binary(&self, a: Matrix2x2, b: Matrix2x2) -> Matrix2x2 {
        Matrix2x2([
            [a.0[0][0] + b.0[0][0], a.0[0][1] + b.0[0][1]],
            [a.0[1][0] + b.0[1][0], a.0[1][1] + b.0[1][1]],
        ])
    }
    
    fn operate(&self, elements: &[Matrix2x2]) -> Matrix2x2 {
        if elements.is_empty() {
            Matrix2x2([[0, 0], [0, 0]]) // zero matrix
        } else {
            elements.iter().fold(Matrix2x2([[0, 0], [0, 0]]), |acc, &x| self.operate_binary(acc, x))
        }
    }
}

impl Group<Matrix2x2> for MatrixAddition {
    fn identity(&self) -> Matrix2x2 {
        Matrix2x2([[0, 0], [0, 0]])
    }
    
    fn inverse(&self, a: Matrix2x2) -> Matrix2x2 {
        Matrix2x2([
            [-a.0[0][0], -a.0[0][1]],
            [-a.0[1][0], -a.0[1][1]],
        ])
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