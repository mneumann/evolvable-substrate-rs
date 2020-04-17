use nalgebra::base::Vector2;
pub use nalgebra::geometry::Point2;

/// The floating-point type used throughout the code
pub type Num = f32;

#[derive(Debug, Copy, Clone)]
pub struct RegionIndex(pub usize);

/// A quadratic region in 2d-space. It can further be divided into four sub-regions
/// (`children`).
#[derive(Debug)]
pub struct Region {
    /// The coordinates of the center point.
    pub center: Point2<Num>,
    /// The distance from the center to the top/bottom and left/right borders.
    /// The width of the region covered by this Region is 2-times as much.
    pub width: Num,
    /// The depth of this node (number of edges between the root and this node).
    pub depth: usize,
    /// The image of the function at (center.x, center.y).
    pub image: Num,
    /// Subdivisions of this Region + variance
    children: Option<(Num, RegionIndex, RegionIndex, RegionIndex, RegionIndex)>,
}

/// Returns the variance^2.
fn variance_sq(a: Num, b: Num, c: Num, d: Num) -> Num {
    let mean = (a + b + c + d) / 4 as Num;
    ((a - mean) * (a - mean))
        + ((b - mean) * (b - mean))
        + ((c - mean) * (c - mean))
        + ((d - mean) * (d - mean))
}

fn top_left() -> Vector2<Num> {
    Vector2::new(-1.0, -1.0)
}

fn top_right() -> Vector2<Num> {
    Vector2::new(1.0, -1.0)
}

fn bottom_left() -> Vector2<Num> {
    Vector2::new(-1.0, 1.0)
}

fn bottom_right() -> Vector2<Num> {
    Vector2::new(1.0, 1.0)
}

#[derive(Debug)]
struct Regions {
    list: Vec<Region>,
}

impl Regions {
    fn new() -> Self {
        Self { list: vec![] }
    }

    fn push(&mut self, region: Region) -> RegionIndex {
        let idx = RegionIndex(self.list.len());
        self.list.push(region);
        idx
    }

    fn at(&self, idx: RegionIndex) -> &Region {
        &self.list[idx.0]
    }

    fn at_mut(&mut self, idx: RegionIndex) -> &mut Region {
        &mut self.list[idx.0]
    }
}

/// Configuration that tells `division_and_initialization` when to stop subdividing a Region.
pub struct PartitionConfig {
    pub min_depth: usize,
    pub max_depth: usize,
    pub variance_sq_threshold: Num,
}

#[derive(Debug)]
pub struct PartitionedSpace {
    regions: Regions,
    root_region_idx: RegionIndex,
}

impl PartitionedSpace {
    // * The papers uses three parameters `outgoing`, `a` and `b`. We simplify this into
    //   the partially-applied function `f`:
    //
    //       f(x, y) = outgoing ? cppn(a, b, x, y) : cppn(x, y, a, b)
    //
    pub fn divide_and_initialize(
        f: &impl Fn(Point2<Num>) -> Num,
        config: &PartitionConfig,
    ) -> Self {
        assert!(config.min_depth <= config.max_depth);
        let mut regions = Regions::new();
        let mut stack = Vec::new();

        let center = Point2::new(0.0, 0.0);
        let root_region = Region {
            center,
            width: 1.0,
            depth: 0,
            image: f(center),
            children: None,
        };

        let root_region_idx = regions.push(root_region);
        stack.push(root_region_idx);

        // contains all the regions that we have to further subdivide

        while let Some(idx) = stack.pop() {
            // divide region with index `idx`
            assert!(regions.at(idx).children.is_none());

            let depth = regions.at(idx).depth;

            let (c1, c2, c3, c4) = regions.at(idx).create_children(f);

            let variance_sq = variance_sq(c1.image, c2.image, c3.image, c4.image);

            let i1 = regions.push(c1);
            let i2 = regions.push(c2);
            let i3 = regions.push(c3);
            let i4 = regions.push(c4);

            regions.at_mut(idx).children = Some((variance_sq, i1, i2, i3, i4));

            if depth < config.max_depth
                && (depth < config.min_depth || variance_sq > config.variance_sq_threshold)
            {
                // further sub-divide
                stack.push(i1);
                stack.push(i2);
                stack.push(i3);
                stack.push(i4);
            }
        }

        Self {
            regions,
            root_region_idx,
        }
    }

    /// Function `select_band_pruned_candidates` has the same basic functionality as
    /// `PruningAndExtraction` from the ES-HyperNEAT paper except that we do not construct
    /// connections here (this can be done within the `select_callback`.
    ///
    /// We want to select those nodes that are enclosed within a "band". For this purpose we either
    /// test leaf nodes or nodes with a low enough variance (< variance_sq_threshold) whether they
    /// are enclosed within a band. If this condition is met, we call the `select_callback`.
    // XXX: Return iterator
    pub fn select_band_pruned_candidates(
        &self,
        f: &impl Fn(Point2<Num>) -> Num,
        variance_sq_threshold: Num,
        band_threshold: Num,
    ) -> Vec<&Region> {
        let mut stack = vec![self.root_region_idx];

        let mut candidates = vec![];

        while let Some(idx) = stack.pop() {
            let region = self.regions.at(idx);
            match region.children {
                None => {
                    // we look at all points that are either a leaf or where the variance is below a threshold (low variance).
                    // for those points, we select those that are located within a "band".
                    if region.is_in_band(f, band_threshold) {
                        candidates.push(region);
                    }
                }
                Some((variance_sq, i1, i2, i3, i4)) => {
                    if variance_sq < variance_sq_threshold {
                        // we look at all points that are either a leaf or where the variance is below a threshold (low variance).
                        // for those points, we select those that are located within a "band".
                        if region.is_in_band(f, band_threshold) {
                            candidates.push(region);
                        }
                    }

                    stack.push(i1);
                    stack.push(i2);
                    stack.push(i3);
                    stack.push(i4);
                }
            }
        }
        candidates
    }

    // XXX: Return iterator
    pub fn each_leaf(&self) -> Vec<&Region> {
        let mut leaves = vec![];
        let mut stack = vec![self.root_region_idx];

        while let Some(idx) = stack.pop() {
            let region = self.regions.at(idx);
            match region.children {
                None => {
                    leaves.push(region);
                }
                Some((_, i1, i2, i3, i4)) => {
                    stack.push(i1);
                    stack.push(i2);
                    stack.push(i3);
                    stack.push(i4);
                }
            }
        }
        leaves
    }
}

impl Region {
    /// Creates a child Region for the given subregion
    fn create_child(&self, f: &impl Fn(Point2<Num>) -> Num, sub_region: Vector2<Num>) -> Self {
        let half_width = self.width / 2.0;
        let center = self.center + (sub_region * half_width);
        Self {
            center,
            width: half_width,
            depth: self.depth + 1,
            image: f(center),
            children: None,
        }
    }

    /// Creates four children, one for each sub region.
    fn create_children(&self, f: &impl Fn(Point2<Num>) -> Num) -> (Region, Region, Region, Region) {
        (
            self.create_child(f, top_left()),
            self.create_child(f, top_right()),
            self.create_child(f, bottom_left()),
            self.create_child(f, bottom_right()),
        )
    }

    fn is_in_horizontal_band(
        &self,
        width: Num,
        f: &impl Fn(Point2<Num>) -> Num,
        band_threshold: Num,
    ) -> bool {
        let image_left = f(self.center + Vector2::new(-width, 0.0));
        let image_right = f(self.center + Vector2::new(width, 0.0));

        let d_left = (self.image - image_left).abs();
        let d_right = (self.image - image_right).abs();

        d_left >= band_threshold && d_right >= band_threshold
    }

    fn is_in_vertical_band(
        &self,
        width: Num,
        f: &impl Fn(Point2<Num>) -> Num,
        band_threshold: Num,
    ) -> bool {
        let image_top = f(self.center + Vector2::new(0.0, -width));
        let image_bottom = f(self.center + Vector2::new(0.0, width));

        let d_top = (self.image - image_top).abs();
        let d_bottom = (self.image - image_bottom).abs();

        d_top >= band_threshold && d_bottom >= band_threshold
    }

    fn is_in_band(&self, f: &impl Fn(Point2<Num>) -> Num, band_threshold: Num) -> bool {
        let width_of_square = self.width * 2.0; // the square has 2-times the width.
        self.is_in_horizontal_band(width_of_square, f, band_threshold)
            || self.is_in_vertical_band(width_of_square, f, band_threshold)
    }
}
