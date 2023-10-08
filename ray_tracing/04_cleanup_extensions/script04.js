let c = document.getElementById("myCanvas");
let ctx = c.getContext("2d");
let width = c.width;
let height = c.height;
var canvasBuffer = ctx.getImageData(0, 0, width, height);
var canvasPitch = canvasBuffer.width * 4;

let vWidth = 1; //viewport width
let vHeight = 1; //viewport height
let projectionPlaneD = 1; // distance between camera and viewport
let bgColor = [0, 0, 0];



let camera = [0, 0, 0];


let sphere1 = new Sphere([0, -1, 3], 1, [255, 0, 0], 500, 0.2); //red
let sphere2 = new Sphere([2, 0, 4], 1, [0, 0, 255], 500, 0.3); // blue
let sphere3 = new Sphere([-2, 0, 4], 1, [0, 255, 0], 10, 0.4); // green
let sphere4 = new Sphere([0, -5001, 0], 5000, [255, 255, 0], 1000, 0.5); // yellow
let sphere5 = new Sphere([0, 2, 10], 1, [255, 0, 255], 500, 0.5); // purple

let spheres = [sphere1, sphere2, sphere3, sphere4, sphere5];

let light1 = new Light("ambient", 0.2, null, null);
let light2 = new Light("point", 0.6, [2, 1, 0], null);
let light3 = new Light("directional", 0.8, null, [1, 4, 4]);

let lights = [light1, light2, light3];

let shiny = false;

let recursionDepth = 3;

function takeScreenshot() {
    c.toBlob((blob) => {
        saveBlob(blob, `screencapture-${c.width}x${c.height}.png`);
      });
}

const saveBlob = (function() {
    const a = document.createElement('a');
    document.body.appendChild(a);
    a.style.display = 'none';
    return function saveData(blob, fileName) {
       const url = window.URL.createObjectURL(blob);
       a.href = url;
       a.download = fileName;
       a.click();
    };
}());

function renderScene() {
    
    // set lighing based on slider values
    setShiny();
    lights[0].intensity = document.getElementById("amb").value / 100;
    lights[1].intensity = document.getElementById("point").value / 100;
    lights[2].intensity = document.getElementById("dir").value / 100;

    // for each pixel in the canvas
    for (let xIndex = (-width/2); xIndex < (width/2); xIndex++) {
        for (let yIndex = (-height/2); yIndex < (height/2); yIndex++) {
            // get viewPort coords // this esssentially details the direction that the camera is looking in
            let D = canvasToViewPort(xIndex, yIndex);
            // find color corresponding to the specific location on the canvas based on raytracing
            let color = traceRay(camera, D, 1, Infinity, recursionDepth);
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
                shiny = true;
            }    
    }
}

function Sphere(center, radius, color, specular, reflective) {
    this.center = center;
    this.radius = radius;
    this.color = color;
    this.specular = specular;
    this.reflective = reflective;
}

function Light(type, intensity, point, direction) {
    this.type = type;
    this.intensity = intensity;
    this.point = point;
    this.direction = direction;
}

// sets the pixel at x,y to color
function putPixel(x, y, color) {
    // map to canvas coordinates
    x = width/2 + x;
    y = height/2 - y - 1;
  
    if (x < 0 || x >= width || y < 0 || y >= height) {
      return;
    }
  
    // not sure what canvas pitch is
    // this just used to offset the canvas buffer to set the appropiate pixel color
    let offset = 4*x + canvasPitch*y;
    canvasBuffer.data[offset++] = color[0];
    canvasBuffer.data[offset++] = color[1];
    canvasBuffer.data[offset++] = color[2];
    canvasBuffer.data[offset++] = 255; // Alpha = 255 (full opacity)
}

// write the canvas buffer to the canvase
function updateCanvas() {
    ctx.putImageData(canvasBuffer, 0, 0);
}

// convert canvas coordinates to the coordinates of the viewport.
// don't quite understand why this is necessary
function canvasToViewPort(x, y) {
    return [(x * vWidth / width), (y * vHeight / height), projectionPlaneD];
}

// get closest intersection sphere from starting point 'origin' along 'direction' between tMin, tMax
function closestIntersection(origin, direction, tMin, tMax) {
    let closest_t = Infinity;
    let closest_sphere = null;

    for (let i = 0; i < spheres.length; i++) {
        let sphere = spheres[i];
        let intersections = intersectRaySphere(origin, direction, sphere);
        let t1 = intersections[0];
        let t2 = intersections[1];

        if (t1 < closest_t && tMin < t1 && t1 < tMax) {
            closest_t = t1;
            closest_sphere = sphere;
        }
        if (t2 < closest_t && tMin < t2 && t2 < tMax) {
            closest_t = t2;
            closest_sphere = sphere;
        }
    }
    return [closest_sphere, closest_t];
}

// trace the ray from viewport location to the clossest sphere.
// return the color of that sphere (or bg color if no sphere along ray)
// to paint the canvas pixel to the correct color in order to display that sphere.
function traceRay(origin, direction, tMin, tMax, recursionDepth) {
    let res = closestIntersection(origin, direction, tMin, tMax);
    let closest_sphere = res[0];
    let closest_t = res[1];

    if (closest_sphere == null) {
        return bgColor;
    }

    let intersectionPoint = computeIntersection(origin, closest_t, direction);
    let normal = subtractPoints(intersectionPoint, closest_sphere.center);
    let len = vectorLength(normal);

    normal = scalarMultiply(1 / len, normal);

    let specular = closest_sphere.specular;
    if (!shiny) {
        specular = -1;
    }
    let intensity = computeLighting(intersectionPoint, normal, scalarMultiply(-1, direction), specular);
    let local_color = scalarMultiply(intensity, closest_sphere.color);

    // if we hit recursion limit, or object is not reflective, we're done.
    let pointReflectivity = closest_sphere.reflective;
    if (recursionDepth <= 0 || pointReflectivity <= 0) {
        return local_color
    }
    // compute reflected color
    let reflectedRay = reflectRay(scalarMultiply(-1, direction), normal) 
    let reflected_color = traceRay(intersectionPoint, reflectedRay, 0.001, Infinity, recursionDepth - 1);

    let local_color1 = scalarMultiply((1 - pointReflectivity), local_color);
    let reflected_color1 = scalarMultiply(pointReflectivity, reflected_color);
    return [local_color1[0] + reflected_color1[0], local_color1[1] + reflected_color1[1], local_color1[2] + reflected_color1[2]];
}

// the function that finds the intersections between a ray and a sphere
function intersectRaySphere(origin, direction, sphere) {
    let co = subtractPoints(origin, sphere.center);
    let a = dotProduct(direction, direction);
    let b = 2 * dotProduct(co, direction);
    let c = dotProduct(co, co) - sphere.radius * sphere.radius;

    let discriminant = (b * b) - (4 * a * c);
    if (discriminant < 0) {
        return [Infinity, Infinity];
    }
    let t1 = (-b + Math.sqrt(discriminant)) / (2 * a);
    let t2 = (-b - Math.sqrt(discriminant)) / (2 * a);

    return [t1, t2];
}

function computeLighting(point, normal, viewDirection, specularity) {
    let intensity = 0.0;
    for (let i = 0; i < lights.length; i++) {
        let light = lights[i];
        let tMax = 1;
        if (light.type == "ambient") {
            intensity += light.intensity;
        } else {
            let lightDirection = [0, 0, 0];
            if (light.type == "point") {
                lightDirection = subtractPoints(light.point, point);
                tMax = 1;
            } else {
                lightDirection = light.direction;
                tMax = Infinity;
            }
            // Shadow check
            let shadow_sphere = closestIntersection(point, lightDirection, 0.001, tMax)[0];
            if (shadow_sphere != null) {
                continue;
            }
            // Diffuse lighting
            let n_dot_l = dotProduct(normal, lightDirection);
            if (n_dot_l > 0) {
                intensity += light.intensity * n_dot_l / (vectorLength(normal) * vectorLength(lightDirection));
            }
            //Specular lighting
            if (specularity != -1) {

                //(2N * <N,L>) - L
                let R = subtractPoints(scalarMultiply(2 * dotProduct(normal, lightDirection), normal), lightDirection);
                let r_dot_v = dotProduct(R, viewDirection);
                if (r_dot_v > 0) {
                    //Li * (<R,V> / |R| * |L|)^ s
                    intensity += light.intensity * Math.pow(r_dot_v / (vectorLength(R) * vectorLength(viewDirection)), specularity);
                }
            }
        }
    }
    return intensity;  
}

function dotProduct(p1, p2) {
    return (p1[0] * p2[0]) + (p1[1] * p2[1]) + (p1[2] * p2[2]);
}

function vectorLength(v) {
    return Math.sqrt(dotProduct(v, v));
}

function scalarMultiply(scalar, point) {
    return [point[0] * scalar, point[1] * scalar, point[2] * scalar];
}

function computeIntersection(origin, t, direction) {
    let tD = scalarMultiply(t, direction);
    return addPoints(origin, tD);
}

function subtractPoints(v1, v2) {
    return [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]];
}

function addPoints(v1, v2) {
    return [v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]];
}

function reflectRay(r, n) {
    let n_dot_r = dotProduct(n, r);
    return subtractPoints(scalarMultiply(2 * n_dot_r, n), r);
}