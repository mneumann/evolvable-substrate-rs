const DELTAS = [
  { x: -1.0, y: -1.0 }, // top left
  { x: 1.0, y: -1.0 }, // top right
  { x: -1.0, y: 1.0 }, // bottom left
  { x: 1.0, y: 1.0 }, // bottom right
];

function mean(regions) {
  let v = 0.0;
  regions.forEach((region) => (v += region.image));
  return v / regions.length;
}

function variance(regions) {
  let _mean = mean(regions);
  let v = 0.0;
  regions.forEach((region) => {
    let dist = region.image - _mean;
    v += dist * dist;
  });
  return Math.sqrt(v);
}

function divide_region(region, config) {
  if (region.children !== null) throw "Assertion";

  let children = config.deltas.map((delta) =>
    create_child_region(region, delta, config.f)
  );
  // v determines how much the image values vary across the child regions.
  // if there is a high variance, there is high information and we want to further sub divide.
  region.variance = variance(children);
  //console.log(v);

  if (region.depth < config.max_depth) {
    if (region.depth < config.min_depth) {
      region.children = children;
    }
    if (region.variance >= config.variance_threshold) {
      region.children = children;
    }
  }

  return region.children;
}

function create_child_region(parent_region, delta, f) {
  let half_width = parent_region.width / 2.0;
  let center = {
    x: parent_region.center.x + delta.x * half_width,
    y: parent_region.center.y + delta.y * half_width,
  };
  return {
    type: "Region",
    center,
    width: half_width,
    depth: parent_region.depth + 1,
    image: f(center),
    children: null,
  };
}

function process_stack(stack, config) {
  while (true) {
    let region = stack.pop();
    if (!region) {
      break;
    }
    let sub_divide = divide_region(region, config);
    if (sub_divide) stack.push(...sub_divide);
  }
}

function eachLeafRegion(region, f) {
  if (!region.children) {
    f(region);
  } else {
    region.children.forEach((child) => eachLeafRegion(child, f));
  }
}

function divide_and_initialize(config) {
  if (config.deltas === undefined) config.deltas = DELTAS;

  if (config.min_depth > config.max_depth) throw "min_depth > max_depth";

  let center = { x: 0, y: 0 };

  let root_region = {
    type: "Region",
    center,
    // this is actual the width from the center to the outside (the width of the square is 2*width)
    width: 1.0,
    depth: 0,
    image: config.f(center),
    children: null,
  };

  process_stack([root_region], config);

  return root_region;
}

window.DivideAndInitialize = divide_and_initialize;

function add_scaled_point(p1, p2, scale) {
  return {
    x: p1.x + scale * p2.x,
    y: p1.y + scale * p2.y,
  };
}

function is_in_linear_band(region, vector, config) {
  let image_left = config.f(add_scaled_point(region.center, vector, -1.0));
  let image_right = config.f(add_scaled_point(region.center, vector, 1.0));
  let d_left = Math.abs(region.image - image_left);
  let d_right = Math.abs(region.image - image_right);

  return d_left >= config.band_threshold && d_right >= config.band_threshold;
}

function is_in_horizontal_band(region, width, config) {
  return is_in_linear_band(region, { x: width, y: 0 }, config);
}

function is_in_vertical_band(region, width, config) {
  return is_in_linear_band(region, { x: 0, y: width }, config);
}

function is_in_band(region, config) {
  let width_of_square = region.width * 2.0; // the square has 2-times the width.
  return (
    is_in_horizontal_band(region, width_of_square, config) ||
    is_in_vertical_band(region, width_of_square, config)
  );
}

/// Function `select_band_pruned_candidates` has the same basic functionality as
/// `PruningAndExtraction` from the ES-HyperNEAT paper except that we do not construct
/// connections here (this can be done within the `select_callback`.
///
/// We want to select those nodes that are enclosed within a "band". For this purpose we either
/// test leaf nodes or nodes with a low enough variance (< variance_sq_threshold) whether they
/// are enclosed within a band. If this condition is met, we call the `select_callback`.
function select_band_pruned_candidates(root_region, config, select_callback) {
  let stack = [root_region];

  while (true) {
    let region = stack.pop();
    if (!region) {
      break;
    }

    // we look at all points that are either a leaf or where the variance is below a threshold (low variance).
    // for those points, we select those that are located within a "band".

    if (region.children) {
      if (
        region.variance < config.variance_threshold &&
        is_in_band(region, config)
      ) {
        select_callback(region);
      }

      stack.push(...region.children);
    } else {
      // leaf
      if (is_in_band(region, config)) {
        select_callback(region);
      }
    }
  }
}

window.select_band_pruned_candidates = select_band_pruned_candidates;

// ------------------------------------
// Plotting support
// ------------------------------------

// ranges are inclusive [min, max]
function make_transformation({ input_domain, output_domain }) {
  let input_width = Math.abs(input_domain.max - input_domain.min);
  let output_width = Math.abs(output_domain.max - output_domain.min);
  let scale = input_width == 0.0 ? 0.0 : output_width / input_width;
  let offset = output_domain.min - input_domain.min * scale;
  return (point) => point * scale + offset;
}

function functions() {
  function zero({ x, y }) {
    return 0.0;
  }

  function circle({ x, y }) {
    let r = x * x + y * y;
    if (r > 0.88) {
      return 0.0;
    } else if (r > 0.22) {
      return 1.0;
    } else if (r > 0.05) {
      return 0.5;
    } else if (r > 0.03) {
      return 0.0;
    } else if (r > 0.01) {
      return 1.0;
    } else {
      return 0.0;
    }
  }

  function four_squares_diag({ x, y }) {
    if ((x < 0.0 && y < 0.0) || (x > 0.0 && y > 0.0)) {
      return 0.0;
    } else {
      return 1.0;
    }
  }

  function sin_cos_pattern({ x, y }) {
    return Math.sin(20.0 * x) * Math.sin(20.0 * y);
  }

  function gauss(u, x) {
    return Math.exp(-Math.pow(x - u, 2) / 4.0) / Math.sqrt(Math.PI * 2.0);
  }

  function s(x) {
    if (Math.abs(x) > 1.0) {
      return 0.0;
    } else {
      return 1.0 - Math.abs(x);
    }
  }

  function figure_c({ x, y }) {
    return gauss(y * 10.0, Math.sin(10.0 * x));
  }

  function figure_e({ x, y }) {
    return (
      s(gauss(10.0 * y, Math.sin(10.0 * x))) +
      s(gauss(0.0, 10.0 * x) * Math.sin(10.0 * y))
    );
  }

  return {
    zero,
    circle,
    four_squares_diag,
    sin_cos_pattern,
    figure_c,
    figure_e,
  };
}

// Math.min(...splat) exceeds stack on Chrome
function Array_min(ary) {
  return ary.reduce((a, b) => Math.min(a, b), ary[0]);
}

function Array_max(ary) {
  return ary.reduce((a, b) => Math.max(a, b), ary[0]);
}

function plot_function({ raster_width, raster_height, f, xrange, yrange }) {
  let pixels = new Float32Array(raster_width * raster_height);

  // raster to image
  let transform_x = make_transformation({
    input_domain: { min: 0, max: raster_width - 1 },
    output_domain: xrange,
  });
  let transform_y = make_transformation({
    input_domain: { min: 0, max: raster_height - 1 },
    output_domain: yrange,
  });

  for (var x = 0; x < raster_width; x++) {
    for (var y = 0; y < raster_height; y++) {
      let tx = transform_x(x);
      let ty = transform_y(y);
      let z = f({ x: tx, y: ty });
      pixels[x + y * raster_width] = z;
    }
  }

  let min = Array_min(pixels);
  let max = Array_max(pixels);

  let normalize_zrange = make_transformation({
    input_domain: { min, max },
    output_domain: { min: 0.0, max: 1.0 },
  });

  let image_to_raster_x = make_transformation({
    output_domain: { min: 0, max: raster_width - 1 },
    input_domain: xrange,
  });
  let image_to_raster_y = make_transformation({
    output_domain: { min: 0, max: raster_height - 1 },
    input_domain: yrange,
  });
  let image_to_raster = function ({ x, y }) {
    return { x: image_to_raster_x(x), y: image_to_raster_y(y) };
  };

  return {
    pixels: pixels.map((z) => normalize_zrange(z)),
    raster_width,
    raster_height,
    image_to_raster,
  };
}

function toImageData({ pixels, raster_width, raster_height }) {
  let rgba_array = new Uint8ClampedArray(4 * raster_width * raster_height);

  for (var i = 0; i < raster_width * raster_height; i++) {
    let pixel = pixels[i];
    rgba_array[i * 4 + 0] = 255 * pixel;
    rgba_array[i * 4 + 1] = 255 * pixel;
    rgba_array[i * 4 + 2] = 255 * pixel;
    rgba_array[i * 4 + 3] = 255;
  }

  return new ImageData(rgba_array, raster_width, raster_height);
}

async function drawCanvas(canvas, config) {
  var ctx = canvas.getContext("2d");
  ctx.save();

  let f = config.f;

  let plot = plot_function({
    raster_width: 400,
    raster_height: 400,
    f,
    xrange: { min: -1.0, max: 1.0 },
    yrange: { min: -1.0, max: 1.0 },
  });

  if (config.show_function) {
    let image_data = toImageData(plot);
    let bmp = await createImageBitmap(image_data);
    ctx.drawImage(bmp, 0, 0, canvas.width, canvas.height);
  } else {
    ctx.fillStyle = "darkgray";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
  }

  ctx.scale(2.0, 2.0);
  let root_region = DivideAndInitialize(config);

  let draw_node = function (node, color_center, color_region) {
    if (node.depth < config.show_from_depth) return;

    let p1 = plot.image_to_raster({
      x: node.center.x - node.width,
      y: node.center.y - node.width,
    });
    let p2 = plot.image_to_raster({
      x: node.center.x + node.width,
      y: node.center.y + node.width,
    });
    let width = p2.x - p1.x;
    let height = p2.y - p1.y;
    let center = plot.image_to_raster(node.center);
    let radius = 0.5 * 0.25 * (12 - node.depth);

    if (config.show_regions) {
      ctx.beginPath();
      ctx.rect(p1.x, p1.y, width, height);
      ctx.lineWidth = 0.25;
      ctx.strokeStyle = color_region;
      ctx.stroke();
    }

    if (config.show_center) {
      ctx.beginPath();
      ctx.arc(center.x, center.y, radius, 0, 2 * Math.PI, false);
      ctx.fillStyle = color_center;
      ctx.fill();
    }
  };

  if (config.show_center || config.show_regions) {
    if (config.show_band_pruned) {
      select_band_pruned_candidates(root_region, config, (node) =>
        draw_node(node, config.color_center, config.color_region)
      );
    } else {
      eachLeafRegion(root_region, (leaf) =>
        draw_node(leaf, config.color_center, config.color_region)
      );
    }
  }
  ctx.restore();
}

window.drawCanvas = drawCanvas;
window.functions = functions;
