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
        a + b
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
        a * b
    }
    
    fn operate(&self, elements: &[i32]) -> i32 {
        if elements.is_empty() {
            1 // identity for multiplication
        } else {
            elements.iter().fold(1, |acc, &x| self.operate_binary(acc, x))
        }
    }
}

// Modular arithmetic (mod 7) group
pub struct Mod7;

impl Algebra<i32> for Mod7 {
    fn operate_binary(&self, a: i32, b: i32) -> i32 {
        (a + b) % 7
    }
    
    fn operate(&self, elements: &[i32]) -> i32 {
        if elements.is_empty() {
            0 // identity for addition mod 7
        } else {
            elements.iter().fold(0, |acc, &x| self.operate_binary(acc, x))
        }
    }
}

impl Group<i32> for Mod7 {
    fn identity(&self) -> i32 {
        0
    }
    
    fn inverse(&self, a: i32) -> i32 {
        (7 - a) % 7
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
