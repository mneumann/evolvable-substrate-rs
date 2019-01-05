/// The floating-point type used throughout the code
type Num = f32;

/// Configuration that tells `division_and_initialization` when to further stop subdividing a QuadPoint.
pub struct Config {
    min_depth: usize,
    max_depth: usize,
    variance_sq_threshold: Num,
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
    /// Returns the variance^2.
    fn variance_sq(children: &[QuadPoint]) -> Num {
        assert!(!children.is_empty());
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
        variance
    }

    pub fn is_leaf(&self) -> bool {
        self.children.is_none()
    }

    pub fn each_leaf(&self, f: &mut impl FnMut(&QuadPoint)) {
        match self.children {
            Some(ref children) => {
                for child in children.iter() {
                    child.each_leaf(f);
                }
            }
            None => {
                f(self);
            }
        }
    }
    pub fn each_node(&self, f: &mut impl FnMut(&QuadPoint)) {
        f(self);
        match self.children {
            Some(ref children) => {
                for child in children.iter() {
                    child.each_node(f);
                }
            }
            None => {}
        }
    }
}

fn create_child_quad_point(
    node: &QuadPoint,
    f: &impl Fn((Num, Num)) -> Num,
    x_mul: Num,
    y_mul: Num,
) -> QuadPoint {
    let x = node.center.x + (x_mul * node.width);
    let y = node.center.y + (y_mul * node.width);
    QuadPoint {
        center: Point2d { x, y },
        width: node.width / 2.0,
        level: node.level + 1,
        f_value: f((x, y)),
        children: None,
    }
}

fn divide_node(
    f: &impl Fn((Num, Num)) -> Num,
    node: &QuadPoint,
    config: &Config,
) -> Option<Box<[QuadPoint; 4]>> {
    let mut children = Box::new([
        create_child_quad_point(&node, f, -0.5, -0.5),
        create_child_quad_point(&node, f, 0.5, -0.5),
        create_child_quad_point(&node, f, -0.5, 0.5),
        create_child_quad_point(&node, f, 0.5, 0.5),
    ]);

    // XXX: variance is not yet know below.
    if node.level < config.min_depth
        || (node.level < config.max_depth
            && QuadPoint::variance_sq(&children[..]) > config.variance_sq_threshold)
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
pub fn division_and_initialization(f: &impl Fn((Num, Num)) -> Num, config: &Config) -> QuadPoint {
    let x = 0.0;
    let y = 0.0;
    let mut root = QuadPoint {
        center: Point2d { x: 0.0, y: 0.0 },
        width: 1.0,
        level: 1,
        f_value: f((x, y)),
        children: None,
    };
    root.children = divide_node(f, &root, config);
    root
}

#[cfg(test)]
mod tests {
    use crate::{division_and_initialization, Config, Num};

    fn f((_x, _y): (Num, Num)) -> Num {
        0.0
    }

    // A similar circular function as given in the ES paper
    fn circle((x, y): (Num, Num)) -> Num {
        let r = x * x + y * y;
        if r > 0.88 {
            0.0
        } else if r > 0.22 {
            1.0
        } else if r > 0.05 {
            0.5
        } else if r > 0.03 {
            0.0
        } else if r > 0.01 {
            1.0
        } else {
            0.0
        }
    }

    fn another((x, y): (Num, Num)) -> Num {
        (20.0 * x).sin() * (20.0 * y).sin()
    }

    fn gauss(u: Num, x: Num) -> Num {
        (-(x - u).powi(2) / 4.0).exp() / (2.0 as f32 * 3.1415).sqrt()
    }

    fn s(x: Num) -> Num {
        if x.abs() > 1.0 {
            0.0
        } else {
            1.0 - x.abs()
        }
    }

    fn figure_c((x, y): (Num, Num)) -> Num {
        gauss(y * 10.0, (10.0 * x).sin())
    }

    fn figure_e((x, y): (Num, Num)) -> Num {
        s(gauss(10.0 * y, (10.0 * x).sin())) + s(gauss(0.0, 10.0 * x) * ((10.0 * y).sin()))
    }

    fn analyse(
        f: impl Fn((f32, f32)) -> f32,
        f2: impl Fn((f32, f32)) -> f32,
        name: &str,
        max_depth: usize,
    ) {
        use splot::{Canvas, ColorPalette, Surface2};
        let surface = Surface2::new(f);
        let gray = ColorPalette::grayscale(256);
        let mut canvas = Canvas::new(1024, 1024);
        canvas.splot(&surface, &gray);
        //let mut img_buffer = surface.plot(&gray, raster);

        let config = Config {
            min_depth: 2,
            max_depth,
            variance_sq_threshold: 0.005,
        };
        let root = division_and_initialization(&f2, &config);
        let transformation = canvas.make_transformation(&(-1.0..=1.0), &(-1.0..=1.0));

        root.each_leaf(&mut |qp| {
            let center = transformation.image_to_raster(qp.center.x, qp.center.y);
            let x1 = qp.center.x - qp.width /* / 2.0 */;
            let y1 = qp.center.y - qp.width /* / 2.0 */;
            let x2 = qp.center.x + qp.width;
            let y2 = qp.center.y + qp.width;
            let r1 = transformation.image_to_raster(x1, y1);
            let r2 = transformation.image_to_raster(x2, y2);
            let w = r2.0 - r1.0;
            canvas.draw_filled_circle(
                center.0,
                center.1,
                (max_depth + 2) as u32 - qp.level as u32,
                (128, 128, 128).into(),
            );
            canvas.draw_square(r1.0, r1.1, w, (255, 0, 0).into());
        });
        canvas.image().save(name).unwrap();
    }

    #[test]
    fn test_circle() {
        analyse(circle, circle, "circle.png", 8);
    }

    #[test]
    fn test_another() {
        analyse(another, another, "another.png", 6);
    }

    #[test]
    fn test_figure_c() {
        analyse(figure_c, figure_c, "figure_c.png", 8);
    }

    #[test]
    fn test_figure_e() {
        analyse(figure_e, figure_e, "figure_e.png", 8);
    }

    #[test]
    fn division_and_initialization_should_stop_at_max_depth() {
        let config = Config {
            min_depth: 1,
            max_depth: 1,
            variance_sq_threshold: 1.0,
        };
        let root = division_and_initialization(&f, &config);
        assert!(!root.is_leaf());
        for i in 0..4 {
            assert!(root.children.as_ref().unwrap()[i].is_leaf());
        }
    }
}
