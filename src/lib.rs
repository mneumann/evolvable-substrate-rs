/// The floating-point type used throughout the code
type Num = f32;

/// Configuration that tells `division_and_initialization` when to further stop subdividing a QuadPoint.
pub struct Config {
    min_depth: usize,
    max_depth: usize,
    variance_threshold: Num,
}

#[derive(Debug, Copy, Clone)]
pub struct Point2d<T> {
    x: T,
    y: T,
}

/// A QuadPoint is a quadratic region in 2d-space which can further be subdivided into four regions (`children`).
#[derive(Debug)]
pub struct QuadPoint {
    /// The coordinates of the center point.
    center: Point2d<Num>,
    /// The width of the region covered by this QuadPoint.
    width: Num,
    /// The depth from root to this node (XXX: actually depth+1).
    level: usize,
    /// The value of f at (center.x, center.y).
    f_value: Num,
    /// Optional subdivisions of this QuadPoint.
    children: Option<Box<[QuadPoint; 4]>>,
}

impl QuadPoint {
    fn variance(&self) -> Num {
        if let Some(ref children) = self.children {
            let mut mean = 0.0;
            for child in children.iter() {
                mean += child.f_value;
            }
            mean /= children.len() as Num;
            let mut variance = 0.0;
            for child in children.iter() {
                let d = child.f_value - mean;
                variance += d * d;
            }
            return variance.sqrt();
        } else {
            0.0
        }
    }

    pub fn is_leaf(&self) -> bool {
        self.children.is_none()
    }
}

fn create_child_quad_point(
    node: &QuadPoint,
    f: &impl Fn(Num, Num) -> Num,
    x_mul: Num,
    y_mul: Num,
) -> QuadPoint {
    let x = node.center.x + (x_mul * node.width);
    let y = node.center.y + (y_mul * node.width);
    QuadPoint {
        center: Point2d { x, y },
        width: node.width / 2.0,
        level: node.level + 1,
        f_value: f(x, y),
        children: None,
    }
}

fn divide_node(
    f: &impl Fn(Num, Num) -> Num,
    node: &QuadPoint,
    config: &Config,
) -> Option<Box<[QuadPoint; 4]>> {
    let mut children = Box::new([
        create_child_quad_point(&node, f, -0.5, -0.5),
        create_child_quad_point(&node, f, 0.5, -0.5),
        create_child_quad_point(&node, f, -0.5, 0.5),
        create_child_quad_point(&node, f, 0.5, 0.5),
    ]);

    if node.level < config.min_depth
        || (node.level < config.max_depth && node.variance() > config.variance_threshold)
    {
        for child in children.iter_mut() {
            child.children = divide_node(f, child, config);
        }
    }

    Some(children)
}

// * The iterative algorithm as presented in the paper does not easily work in Rust
//   due to aliasing. We use a recursive approach instead. A recursive approach is
//   ok, as recusion depth is usually very low (< 20).
//
// * The papers uses three parameters `outgoing`, `a` and `b`. We simplify this into
//   the partially-applied function `f`:
//
//       f(x, y) = outgoing ? cppn(a, b, x, y) : cppn(x, y, a, b)
//
pub fn division_and_initialization(f: &impl Fn(Num, Num) -> Num, config: &Config) -> QuadPoint {
    let x = 0.0;
    let y = 0.0;
    let mut root = QuadPoint {
        center: Point2d { x, y },
        width: 1.0,
        level: 1,
        f_value: f(x, y),
        children: None,
    };
    root.children = divide_node(f, &root, config);
    root
}

#[cfg(test)]
mod tests {
    use crate::{division_and_initialization, Config, Num};

    fn f(_x: Num, _y: Num) -> Num {
        return 0.0;
    }

    #[test]
    fn division_and_initialization_should_stop_at_max_depth() {
        let config = Config {
            min_depth: 1,
            max_depth: 1,
            variance_threshold: 1.0,
        };
        let root = division_and_initialization(&f, &config);
        assert!(!root.is_leaf());
        for i in 0..4 {
            assert!(root.children.as_ref().unwrap()[i].is_leaf());
        }
    }
}
