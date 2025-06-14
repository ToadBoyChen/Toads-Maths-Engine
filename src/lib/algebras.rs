// mod libs {
//     include!("../lib/specialised_algebras.rs");
// }

// use crate::libs::*;

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