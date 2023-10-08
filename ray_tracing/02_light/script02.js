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


let sphere1 = new Sphere(new Point(0, -1, 3), 1, new Color(255, 0, 0), 500); //red
let sphere2 = new Sphere(new Point(2, 0, 4), 1, new Color(0, 0, 255), 500); // blue
let sphere3 = new Sphere(new Point(-2, 0, 4), 1, new Color(0, 255, 0), 10); // green
let sphere4 = new Sphere(new Point(0, -5001, 0), 5000, new Color(255, 255, 0), 1000); // yellow
let sphere5 = new Sphere(new Point(0, 2, 10), 1, new Color(255, 0, 255), 500); // purple

let spheres = [sphere1, sphere2, sphere3, sphere4, sphere5];

let light1 = new Light("ambient", 0.2, null, null);
let light2 = new Light("point", 0.6, new Point(2, 1, 0), null);
let light3 = new Light("directional", 0.8, null, new Point(1, 4, 4));

let lights = [light1, light2, light3];

let shiny = false;


function renderScene() {

    // set lighing based on slider values
    setShiny();
    lights[0].intensity = document.getElementById("amb").value / 100;
    lights[1].intensity = document.getElementById("point").value / 100;
    lights[2].intensity = document.getElementById("dir").value / 100;

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
}

renderScene();

function setShiny() {
    var ele = document.getElementsByName("surface");

    for (i = 0; i < ele.length; i++) {
        if (ele[i].checked)
            if (ele[i].value == "shiny") {
                shiny = true;
            }
            else {
                shiny = false;
            }
    }
}



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

function Sphere(center, radius, color, specular) {
    this.center = center;
    this.radius = radius;
    this.color = color;
    this.specular = specular;
}

function Light(type, intensity, point, direction) {
    this.type = type;
    this.intensity = intensity;
    this.point = point;
    this.direction = direction;
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
    return new Point((x * v_width / width), (y * v_height / height), projection_plane_d);
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

    let P = compute_intersection(cam, closest_t, D);
    let N = subtract_points(P, closest_sphere.center);
    let len = vector_length(N);
    N.x = N.x / len;
    N.y = N.y / len;
    N.z = N.z / len;

    let specular = closest_sphere.specular;
    if (!shiny) {
        specular = -1;
    }
    let intensity = computeLighting(P, N, scalar_multiply(-1, D), specular);

    return apply_color_intensity(closest_sphere.color, intensity);
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

function computeLighting(P, N, V, s) {
    let intensity = 0.0;
    for (let i = 0; i < lights.length; i++) {
        let light = lights[i];
        if (light.type == "ambient") {
            intensity += light.intensity;
        } else {
            let L = new Point(0, 0, 0);
            if (light.type == "point") {
                L.x = light.point.x - P.x;
                L.y = light.point.y - P.y;
                L.z = light.point.z - P.z;
            } else {
                L = light.direction;
            }
            // Diffuse lighting
            let n_dot_l = dot_product(N, L);
            if (n_dot_l > 0) {
                intensity += light.intensity * n_dot_l / (vector_length(N) * vector_length(L));
            }
            //Specular lighting
            if (s != -1) {

                //(2N * <N,L>) - L
                let R = subtract_points(scalar_multiply(2 * dot_product(N, L), N), L);
                let r_dot_v = dot_product(R, V);
                if (r_dot_v > 0) {
                    //Li * (<R,V> / |R| * |L|)^ s
                    intensity += light.intensity * Math.pow(r_dot_v / (vector_length(R) * vector_length(V)), s);
                }
            }
        }
    }
    return intensity;
}

function dot_product(p1, p2) {
    return (p1.x * p2.x) + (p1.y * p2.y) + (p1.z * p2.z);
}

function vector_length(v) {
    return Math.sqrt(dot_product(v, v));
}

function scalar_multiply(scalar, point) {
    let new_point = new Point(point.x * scalar, point.y * scalar, point.z * scalar);
    return new_point;
}

function compute_intersection(O, t, D) {
    D.x = t * D.x;
    D.y = t * D.y;
    D.z = t * D.z;

    return new Point(O.x + D.x, O.y + D.y, O.z + D.z);
}

function subtract_points(v1, v2) {
    return new Point(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

function apply_color_intensity(color, i) {
    let color2 = new Color(color.r * i, color.g * i, color.b * i);
    return color2;
}
