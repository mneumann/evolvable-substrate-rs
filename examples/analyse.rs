use evolvable_substrate::{Num, PartitionConfig, PartitionedSpace, Point2, Region};

/// Some test image functions
mod fns {
    use super::{Num, Point2};

    // Four squares where diagonal squares have the same color.
    pub fn four_squares_diag(p: Point2<Num>) -> Num {
        if p.x < 0.0 && p.y < 0.0 || p.x > 0.0 && p.y > 0.0 {
            0.0
        } else {
            1.0
        }
    }

    pub fn zero(_: Point2<Num>) -> Num {
        0.0
    }

    // A similar circular function as given in the ES paper
    pub fn circle(p: Point2<Num>) -> Num {
        let r = p.x * p.x + p.y * p.y;
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

    pub fn another(p: Point2<Num>) -> Num {
        (20.0 * p.x).sin() * (20.0 * p.y).sin()
    }

    pub fn gauss(u: Num, x: Num) -> Num {
        (-(x - u).powi(2) / 4.0).exp() / (std::f32::consts::PI * 2.0).sqrt()
    }

    pub fn s(x: Num) -> Num {
        if x.abs() > 1.0 {
            0.0
        } else {
            1.0 - x.abs()
        }
    }

    pub fn figure_c(p: Point2<Num>) -> Num {
        gauss(p.y * 10.0, (10.0 * p.x).sin())
    }

    pub fn figure_e(p: Point2<Num>) -> Num {
        s(gauss(10.0 * p.y, (10.0 * p.x).sin())) + s(gauss(0.0, 10.0 * p.x) * ((10.0 * p.y).sin()))
    }
}

fn make_canvas(f: impl Clone + Fn(Point2<f32>) -> f32) -> splot::Canvas {
    use splot::{Canvas, ColorPalette, Surface2};
    let surface = Surface2::new(|(x, y)| f(Point2::new(x, y)));
    let gray = ColorPalette::grayscale(256);
    let mut canvas = Canvas::new(1024, 1024);
    canvas.splot(&surface, &gray);
    canvas
}

fn draw_region(
    max_depth: usize,
    qp: &Region,
    canvas: &mut splot::Canvas,
    transformation: &splot::canvas::RasterTransformation,
    color: (u8, u8, u8),
) {
    let center = transformation.image_to_raster(qp.center.x, qp.center.y);

    let x1 = qp.center.x - qp.width;
    let y1 = qp.center.y - qp.width;
    let x2 = qp.center.x + qp.width;
    let y2 = qp.center.y + qp.width;
    let r1 = transformation.image_to_raster(x1, y1);
    let r2 = transformation.image_to_raster(x2, y2);
    let w = r2.0 - r1.0;

    let radius = (max_depth + 2) as u32 - qp.depth as u32;
    canvas.draw_filled_circle(center.0, center.1, radius, color.into());
    canvas.draw_square(r1.0, r1.1, w, color.into());
}

fn analyse(
    f: impl Clone + Fn(Point2<f32>) -> f32,
    name: &str,
    min_depth: usize,
    max_depth: usize,
    variance_sq_threshold: f32,
    prune_variance_sq_threshold: f32,
    band_threshold: f32,
) {
    let config = PartitionConfig {
        min_depth,
        max_depth,
        variance_sq_threshold,
    };
    let root = PartitionedSpace::divide_and_initialize(&f, &config);

    let mut canvas = make_canvas(&f);
    let transformation = canvas.make_transformation(&(-1.0..=1.0), &(-1.0..=1.0));
    for region in root.each_leaf() {
        draw_region(max_depth, region, &mut canvas, &transformation, (0, 255, 0));
    }
    for region in
        root.select_band_pruned_candidates(&f, prune_variance_sq_threshold, band_threshold)
    {
        draw_region(max_depth, region, &mut canvas, &transformation, (255, 0, 0));
    }

    canvas.image().save(format!("{}.png", name)).unwrap();
}

fn main() {
    analyse(fns::circle, "circle", 3, 7, 0.2 * 0.2, 0.2 * 0.2, 0.5);

    analyse(fns::circle, "circle2", 4, 8, 0.2 * 0.2, 0.1 * 0.1, 0.9);

    analyse(fns::another, "another", 2, 6, 0.1, 0.5, 0.01);
    analyse(fns::another, "another2", 2, 6, 0.1, 0.5, 0.1);

    analyse(fns::figure_c, "figure_c", 2, 8, 0.1, 0.05, 0.1);

    analyse(fns::figure_e, "figure_e", 2, 8, 0.2 * 0.2, 0.4 * 0.4, 0.1);
    analyse(fns::figure_e, "figure_e2", 2, 8, 0.2 * 0.2, 0.1 * 0.1, 0.04);

    analyse(
        fns::four_squares_diag,
        "four_squares_diag",
        0,
        2,
        0.0,
        0.0,
        0.1,
    );
}
