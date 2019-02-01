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

/// A QuadPoint is a quadratic region in 2d-space. It can further be divided into four sub-regions
/// (`children`).
#[derive(Debug)]
pub struct QuadPoint {
    /// The coordinates of the center point.
    center: Point2d<Num>,
    /// The distance from the center to the top/bottom and left/right borders.
    /// The width of the region covered by this QuadPoint is 2-times as much.
    width: Num,
    /// The depth of this node (number of edges between the root and this node).
    depth: usize,
    /// The image of the function at (center.x, center.y).
    image: Num,
    /// Subdivisions of this QuadPoint.
    children: Vec<QuadPoint>,
}

/// Returns the variance^2.
fn variance_sq(children: &[QuadPoint]) -> Num {
    let mean = children.iter().map(|c| c.image).sum::<Num>() / children.len() as Num;
    children.iter().map(|c| c.image - mean).map(|d| d * d).sum()
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum SubRegion {
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight,
}

impl QuadPoint {
    pub fn children(&self) -> &[QuadPoint] {
        &self.children[..]
    }

    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    /// Returns Infinity if children is empty
    pub fn variance_sq(&self) -> Num {
        variance_sq(self.children())
    }

    /// Creates a child QuadPoint for the given SubRegion
    fn create_child(&self, f: &impl Fn((Num, Num)) -> Num, sub_region: SubRegion) -> Self {
        let (x_dir, y_dir) = match sub_region {
            SubRegion::TopLeft => (-1.0, -1.0),
            SubRegion::TopRight => (1.0, -1.0),
            SubRegion::BottomLeft => (-1.0, 1.0),
            SubRegion::BottomRight => (1.0, 1.0),
        };

        let half_width = self.width / 2.0;
        let x = self.center.x + (x_dir * half_width);
        let y = self.center.y + (y_dir * half_width);
        Self {
            center: Point2d { x, y },
            width: half_width,
            depth: self.depth + 1,
            image: f((x, y)),
            children: Vec::new(),
        }
    }

    /// Creates four children, one for each SubRegion.
    fn create_children(&self, f: &impl Fn((Num, Num)) -> Num) -> Vec<QuadPoint> {
        vec![
            self.create_child(f, SubRegion::TopLeft),
            self.create_child(f, SubRegion::TopRight),
            self.create_child(f, SubRegion::BottomLeft),
            self.create_child(f, SubRegion::BottomRight),
        ]
    }

    // * The iterative algorithm as presented in the paper does not easily work in Rust
    //   due to aliasing. We use a recursive approach instead. The recursive approach works
    //   well as recusion depth is usually very low (< 20).
    //
    // * The papers uses three parameters `outgoing`, `a` and `b`. We simplify this into
    //   the partially-applied function `f`:
    //
    //       f(x, y) = outgoing ? cppn(a, b, x, y) : cppn(x, y, a, b)
    //
    pub fn division_and_initialization(
        f: &impl Fn((Num, Num)) -> Num,
        config: &Config,
    ) -> QuadPoint {
        let x = 0.0;
        let y = 0.0;
        let mut root = QuadPoint {
            center: Point2d { x, y },
            width: 1.0,
            depth: 0,
            image: f((x, y)),
            children: Vec::new(),
        };
        root.divide_rec(f, 0, config);
        root
    }

    fn is_in_horizontal_band(
        &self,
        width: Num,
        f: &impl Fn((Num, Num)) -> Num,
        band_threshold: Num,
    ) -> bool {
        let image_left = f((self.center.x - width, self.center.y));
        let image_right = f((self.center.x + width, self.center.y));

        let d_left = (self.image - image_left).abs();
        let d_right = (self.image - image_right).abs();

        d_left >= band_threshold && d_right >= band_threshold
    }

    fn is_in_vertical_band(
        &self,
        width: Num,
        f: &impl Fn((Num, Num)) -> Num,
        band_threshold: Num,
    ) -> bool {
        let image_top = f((self.center.x, self.center.y - width));
        let image_bottom = f((self.center.x, self.center.y + width));

        let d_top = (self.image - image_top).abs();
        let d_bottom = (self.image - image_bottom).abs();

        d_top >= band_threshold && d_bottom >= band_threshold
    }

    fn is_in_band(&self, f: &impl Fn((Num, Num)) -> Num, band_threshold: Num) -> bool {
        let width_of_square = self.width * 2.0; // the square has 2-times the width.
        self.is_in_horizontal_band(width_of_square, f, band_threshold)
            || self.is_in_vertical_band(width_of_square, f, band_threshold)
    }

    /// Recursively divide this node and child nodes.
    pub(crate) fn divide_rec(
        &mut self,
        f: &impl Fn((Num, Num)) -> Num,
        depth: usize,
        config: &Config,
    ) {
        assert!(self.is_leaf());
        debug_assert!(config.min_depth <= config.max_depth);

        if depth < config.max_depth {
            self.children = self.create_children(f);
            if depth < config.min_depth || self.variance_sq() > config.variance_sq_threshold {
                for child in self.children.iter_mut() {
                    child.divide_rec(f, depth + 1, config);
                }
            }
        }
    }

    /// Function `select_band_pruned_candidates` has the same basic functionality as
    /// `PruningAndExtraction` from the ES-HyperNEAT paper except that we do not construct
    /// connections here (this can be done within the `select_callback`.
    ///
    /// We want to select those nodes that are enclosed within a "band". For this purpose we either
    /// test leaf nodes or nodes with a low enough variance (< variance_sq_threshold) whether they
    /// are enclosed within a band. If this condition is met, we call the `select_callback`.
    pub fn select_band_pruned_candidates(
        &self,
        f: &impl Fn((Num, Num)) -> Num,
        variance_sq_threshold: Num,
        band_threshold: Num,
        select_callback: &mut impl FnMut(&QuadPoint),
    ) {
        for child in self.children.iter() {
            if child.children.is_empty() || child.variance_sq() < variance_sq_threshold {
                // we look at all points that are either a leaf of where the variance is below a threshold (low variance).
                // for those points, we select those that are located within a "band".
                let is_candidate = child.is_in_band(f, band_threshold);
                if is_candidate {
                    select_callback(child);
                }
            } else {
                child.select_band_pruned_candidates(
                    f,
                    variance_sq_threshold,
                    band_threshold,
                    select_callback,
                );
            }
        }
    }

    pub fn each_leaf(&self, f: &mut impl FnMut(&QuadPoint)) {
        if self.children.is_empty() {
            f(self);
        } else {
            for child in self.children.iter() {
                child.each_leaf(f);
            }
        }
    }

    pub fn each_node(&self, f: &mut impl FnMut(&QuadPoint)) {
        f(self);
        for child in self.children.iter() {
            child.each_node(f);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::Config;

    /// Some test image functions
    pub mod fns {
        use crate::Num;

        // Four squares where diagonal squares have the same color.
        pub fn four_squares_diag((x, y): (Num, Num)) -> Num {
            if x < 0.0 && y < 0.0 || x > 0.0 && y > 0.0 {
                0.0
            } else {
                1.0
            }
        }

        pub fn zero((_x, _y): (Num, Num)) -> Num {
            0.0
        }

        // A similar circular function as given in the ES paper
        pub fn circle((x, y): (Num, Num)) -> Num {
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

        pub fn another((x, y): (Num, Num)) -> Num {
            (20.0 * x).sin() * (20.0 * y).sin()
        }

        fn gauss(u: Num, x: Num) -> Num {
            (-(x - u).powi(2) / 4.0).exp() / (std::f32::consts::PI * 2.0).sqrt()
        }

        fn s(x: Num) -> Num {
            if x.abs() > 1.0 {
                0.0
            } else {
                1.0 - x.abs()
            }
        }

        pub fn figure_c((x, y): (Num, Num)) -> Num {
            gauss(y * 10.0, (10.0 * x).sin())
        }

        pub fn figure_e((x, y): (Num, Num)) -> Num {
            s(gauss(10.0 * y, (10.0 * x).sin())) + s(gauss(0.0, 10.0 * x) * ((10.0 * y).sin()))
        }

    }

    #[test]
    fn test_variance_zero() {
        use crate::{Point2d, QuadPoint};

        let f = fns::zero;
        let x = 0.0;
        let y = 0.0;
        let mut root = QuadPoint {
            center: Point2d { x, y },
            width: 1.0,
            depth: 0,
            image: f((x, y)),
            children: Vec::new(),
        };

        // A quad point with no children has no variance (we cannot calculate it!)
        assert!(root.is_leaf());

        // If we divide a quadpoint, the quadpoint will have children
        root.divide_rec(
            &f,
            0,
            &Config {
                min_depth: 1,
                max_depth: 1,
                variance_sq_threshold: 0.0,
            },
        );

        // And as such, it will have a variance.
        assert_eq!(0.0, crate::variance_sq(&root.children()));
    }

    #[test]
    fn test_variance_four_squares_diag() {
        use crate::{Point2d, QuadPoint};

        let f = fns::four_squares_diag;
        let x = 0.0;
        let y = 0.0;
        let mut root = QuadPoint {
            center: Point2d { x, y },
            width: 1.0,
            depth: 0,
            image: f((x, y)),
            children: Vec::new(),
        };

        // A quad point with no children has no variance (we cannot calculate it!)
        assert!(root.is_leaf());

        // If we divice a quadpoint, the quadpoint will have children
        root.divide_rec(
            &f,
            0,
            &Config {
                min_depth: 1,
                max_depth: 1,
                variance_sq_threshold: 0.0,
            },
        );

        // And as such, it will have a variance.
        // mean = 0.5, variance**2 = 4 * (0.5**2) = 4 * 0.25 = 1
        assert_eq!(1.0, crate::variance_sq(&root.children()));
    }

    fn analyse(
        f: impl Clone + Fn((f32, f32)) -> f32,
        name: &str,
        min_depth: usize,
        max_depth: usize,
        variance_sq_threshold: f32,
        prune_variance_sq_threshold: f32,
        band_threshold: f32,
        prune: bool,
    ) {
        use crate::QuadPoint;
        use splot::{Canvas, ColorPalette, Surface2};
        let surface = Surface2::new(f.clone());
        let gray = ColorPalette::grayscale(256);
        let mut canvas = Canvas::new(1024, 1024);
        canvas.splot(&surface, &gray);

        let config = Config {
            min_depth,
            max_depth,
            variance_sq_threshold,
        };
        let root = QuadPoint::division_and_initialization(&f, &config);

        let transformation = canvas.make_transformation(&(-1.0..=1.0), &(-1.0..=1.0));

        let draw = &mut |qp: &QuadPoint| {
            let center = transformation.image_to_raster(qp.center.x, qp.center.y);

            let x1 = qp.center.x - qp.width;
            let y1 = qp.center.y - qp.width;
            let x2 = qp.center.x + qp.width;
            let y2 = qp.center.y + qp.width;
            let r1 = transformation.image_to_raster(x1, y1);
            let r2 = transformation.image_to_raster(x2, y2);
            let w = r2.0 - r1.0;
            let color = (255, 0, 0);

            let radius = (max_depth + 2) as u32 - qp.depth as u32;
            canvas.draw_filled_circle(center.0, center.1, radius, color.into());

            canvas.draw_square(r1.0, r1.1, w, (255, 0, 0).into());
        };

        if prune {
            root.select_band_pruned_candidates(
                &f,
                prune_variance_sq_threshold,
                band_threshold,
                draw,
            );
        } else {
            root.each_leaf(draw);
        }

        canvas.image().save(name).unwrap();
    }

    #[test]
    fn test_circle() {
        analyse(
            fns::circle,
            "circle.png",
            3,
            7,
            0.2 * 0.2,
            0.2 * 0.2,
            0.5,
            false,
        );
        analyse(
            fns::circle,
            "circle_prune.png",
            3,
            7,
            0.2 * 0.2,
            0.1 * 0.1,
            0.9,
            true,
        );
    }

    #[test]
    fn test_another() {
        analyse(fns::another, "another.png", 2, 6, 0.1, 0.5, 0.01, false);
        analyse(fns::another, "another_prune.png", 2, 6, 0.1, 0.5, 0.1, true);
    }

    #[test]
    fn test_figure_c() {
        analyse(fns::figure_c, "figure_c.png", 2, 8, 0.1, 0.05, 0.1, false);
    }

    #[test]
    fn test_figure_e() {
        analyse(
            fns::figure_e,
            "figure_e.png",
            2,
            8,
            0.2 * 0.2,
            0.4 * 0.4,
            0.1,
            false,
        );
        analyse(
            fns::figure_e,
            "figure_e_prune.png",
            2,
            8,
            0.2 * 0.2,
            0.1 * 0.1,
            0.04,
            true,
        );
    }

    #[test]
    fn test_four_squares_diag() {
        analyse(
            fns::four_squares_diag,
            "four_squares_diag.png",
            0,
            2,
            0.0,
            0.0,
            0.1,
            false,
        );
    }

    #[test]
    fn division_and_initialization_should_stop_at_max_depth() {
        use crate::QuadPoint;
        let config = Config {
            min_depth: 1,
            max_depth: 1,
            variance_sq_threshold: 1.0,
        };
        let root = QuadPoint::division_and_initialization(&fns::zero, &config);
        assert!(!root.is_leaf());
        for i in 0..4 {
            assert!(root.children()[i].is_leaf());
        }
    }
}
