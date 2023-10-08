let c = document.getElementById("myCanvas");
let ctx = c.getContext("2d");
let width = c.width;
let height = c.height;
var canvas_buffer = ctx.getImageData(0, 0, width, height);
var canvas_pitch = canvas_buffer.width * 4;

let v_width = 1; //viewport width
let v_height = 1; //viewport height
let projection_plane_d = 1; // distance between camera and viewport
let bg_color = new Color(255, 255, 255);



let camera = new Point(0, 0, 0);

let sphere1 = new Sphere(new Point(0, -1, 3), 1, new Color(255, 0, 0)); //red
let sphere2 = new Sphere(new Point(2, 0, 3), 1, new Color(0, 0, 255)); // blue
let sphere3 = new Sphere(new Point(-2, 0, 4), 1, new Color(0, 255, 0)); // green

let spheres = [sphere1, sphere2, sphere3];

// for each pixel in the canvas
for (let xIndex = (-width/2); xIndex < (width/2); xIndex++) {
    for (let yIndex = (-height/2); yIndex < (height/2); yIndex++) {
        // get viewPort coords
        let D = canvasToViewPort(xIndex, yIndex);
        // find color corresponding to the specific location on the canvas based on raytracing
        let color = traceRay(camera, D, 1, Infinity);
        // place color at that pixel
        putPixel(xIndex, yIndex, color);
    }
}

// draw the canvas buffer
updateCanvas();


function Color(r, g, b) {
    this.r = r;
    this.g = g;
    this.b = b;
}

function Point(x, y, z) {
    this.x = x;
    this.y = y;
    this.z = z;
}

function Sphere(center, radius, color) {
    this.center = center;
    this.radius = radius;
    this.color = color;
}

// the function that writes a pixes to the canvas buffer
function putPixel(x, y, color) {
    // map to canvas coordinates
    x = width/2 + x;
    y = height/2 - y - 1;

    if (x < 0 || x >= width || y < 0 || y >= height) {
      return;
    }

    // not sure what canvas pitch is
    let offset = 4*x + canvas_pitch*y;
    canvas_buffer.data[offset++] = color.r;
    canvas_buffer.data[offset++] = color.g;
    canvas_buffer.data[offset++] = color.b;
    canvas_buffer.data[offset++] = 255; // Alpha = 255 (full opacity)
}

// write the canvas buffer to the canvase
function updateCanvas() {
    ctx.putImageData(canvas_buffer, 0, 0);
}

// convert canvas coordinates to the coordinates of the viewport.
// don't quite understand why this is necessary
function canvasToViewPort(x, y) {
    return new
	((x * v_width / width), (y * v_height / height), projection_plane_d);
}

// trace the ray from viewport location to the clossest sphere.
// return the color of that sphere (or bg color if no sphere along ray)
// to paint the canvas pixel to the correct color in order to display that sphere.
function traceRay(cam, D, t_min, t_max) {
    let closest_t = Infinity;
    let closest_sphere = null;

    for (let i = 0; i < spheres.length; i++) {
        let sphere = spheres[i];
        let intersections = intersectRaySphere(cam, D, sphere);
        let t1 = intersections[0];
        let t2 = intersections[1];


        if (t1 < closest_t && t_min < t1 && t1 < t_max) {
            closest_t = t1;
            closest_sphere = sphere;
        }
        if (t2 < closest_t && t_min < t2 && t2 < t_max) {
            closest_t = t2;
            closest_sphere = sphere;
        }
    }
    if (closest_sphere == null) {
        return bg_color;
    }

    return closest_sphere.color;
}

// the function that finds the intersections between a ray and a sphere
function intersectRaySphere(cam, D, sphere) {
    let rad = sphere.radius;
    let center = sphere.center;
    let co = new Point(cam.x - center.x, cam.y - center.y, cam.z - center.z);
    let a = dot_product(D, D);
    let b = 2 * dot_product(co, D);
    let c = dot_product(co, co) - rad * rad;

    discriminant = (b * b) - (4 * a * c);
    if (discriminant < 0) {
        return [Infinity, Infinity];
    }
    let t1 = (-b + Math.sqrt(discriminant)) / (2 * a);
    let t2 = (-b - Math.sqrt(discriminant)) / (2 * a);

    return [t1, t2];
}

function dot_product(p1, p2) {
    return (p1.x * p2.x) + (p1.y * p2.y) + (p1.z * p2.z);
}
