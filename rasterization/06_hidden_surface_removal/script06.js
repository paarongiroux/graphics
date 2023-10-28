let c = document.getElementById("myCanvas");
let ctx = c.getContext("2d");
let width = c.width;
let height = c.height;
let canvas_buffer = ctx.getImageData(0, 0, width, height);
let canvas_pitch = canvas_buffer.width * 4;

let v_width = 1; //viewport width
let v_height = 1; //viewport height
let projection_plane_d = 1; // distance between camera and viewport

let camera = [0,0,0];
let cameraRotation = [0, 0, 0];

let black = [0, 0, 0];
let green = [0, 255, 0];
let blue = [0, 0, 255];
let red = [255 , 0, 0];
let yellow = [255, 255, 0];
let purple = [255, 0, 255];
let cyan = [0, 255, 255];

let depthBuffer = []; // x is width , y is height so dimensions should be width, height
initializeDepthBuffer();


//TO-DO: debug clipping.
    // seems to be clipping whole objects and not doing the triangle clipping properly
    // top and bottom planes are easy to test with

// Add event listener on keypress
document.addEventListener('keypress', (event) => {
    var name = event.key;
    var code = event.code;
    cameraRotation[1] = cameraRotation[1] % 360;
    // Alert the key name and key code on keydown
    if (code==="KeyW"){
        let mov = computeCamMovement(cameraRotation[1]);
        camera[2] += mov[0];
        camera[0] += mov[1];
        renderScene();
        updateCanvas();
    }
    if (code==="KeyS"){
        let mov = computeCamMovement(cameraRotation[1]);
        camera[2] -= mov[0];
        camera[0] -= mov[1];
        renderScene();
        updateCanvas();
    }
    if (code==="KeyA"){
        let mov = computeCamMovement(cameraRotation[1] - 90);
        camera[2] += mov[0];
        camera[0] += mov[1];
        //angle -90
        // camera[0] -= 1;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyD"){
        let mov = computeCamMovement(cameraRotation[1] + 90);
        camera[2] += mov[0];
        camera[0] += mov[1];
        // angle +90
        // camera[0] += 1;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyQ"){
        camera[1] += 1;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyE"){
        camera[1] -= 1;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyI"){
        cameraRotation[0] += 2;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyO"){
        cameraRotation[0] -= 2;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyK"){
        cameraRotation[1] -= 5;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyL"){
        cameraRotation[1] += 5;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyN"){
        cameraRotation[2] += 2;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyM"){
        cameraRotation[2] -= 2;
        renderScene();
        updateCanvas();
    }
  }, false);


let cubeVertices = [
    [1, 1, 1],
    [-1, 1, 1],
    [-1, -1, 1],
    [1, -1, 1],
    [1, 1, -1],
    [-1, 1, -1],
    [-1, -1, -1],
    [1, -1, -1]
];

let cubeTriangles = [
    [0, 1, 2, red],
    [0, 2, 3, red],
    [4, 0, 3, yellow],
    [4, 3, 7, yellow],
    [5, 4, 7, green],
    [5, 7, 6, green],
    [1, 5, 6, blue],
    [1, 6, 2, blue],
    [4, 5, 1, purple],
    [4, 1, 0, purple],
    [2, 6, 7, cyan],
    [2, 7, 3, cyan]
];

let triangleVertices = [
    [1, 1, 1],
    [-1, 1, 1],
    [-1, -1, 1]
];

let triangleTriangles = [
    [0,1,2, red]
];


let cubeModel = new Model(cubeVertices, cubeTriangles);
let triangleModel = new Model(triangleVertices, triangleTriangles);

let cube1 = new Instance(cubeModel, new Transform(1, [0, 0, 0], [-1.5, 0, 7]));
let cube2 = new Instance(cubeModel, new Transform(1, [45, 45, 45], [1.5, 0, 10]));
let cube3 = new Instance(cubeModel, new Transform(1, [0, 0, 0], [-1.5, 0, -7]));

let triangle1 = new Instance(triangleModel, new Transform(1, [0,0,0], [-1.5, 0, 7]));

let scene = [cube1, cube2, cube3];
let s2 = Math.sqrt(2);
let planeNear = new Plane([0, 0, 1], -projection_plane_d, "NEAR");
let planeLeft = new Plane([s2, 0, s2], 0, "LEFT");
let planeRight = new Plane([-s2, 0, s2], 0, "RIGHT");
let planeTop = new Plane([0, -s2, s2], 0, "TOP");
let planeBottom = new Plane([0, s2, s2], 0, "BOTTOM");

let clippingPlanes = [planeNear, planeLeft, planeRight, planeTop, planeBottom];



renderScene();

updateCanvas();

function Model(vertices, triangles) {
    this.vertices = vertices;
    this.triangles = triangles;
}

function Instance(model, transform) {
    this.model = model;
    this.transform = transform;
}

function Transform(scale, rotation, position) {
    this.scale = scale;
    this.rotation = rotation; // [x,y,z] degrees to rotate about each axis
    this.position = position;
}

function Plane(normal, d, name) {
    this.normal = normal;
    this.d = d;
    this.name = name;
}

function initializeDepthBuffer() {
    for (let i = 0; i < width; i++) {
        let tempArray = Array(height).fill(Infinity);
        depthBuffer[i] = tempArray;
    }
}


function clipInstance(instance, planes) {
    let realTrianglesInstance = new Model(instance.vertices, instance.triangles);
    let triangles = createRealTrianglesList(instance);
    realTrianglesInstance.triangles = triangles;
    for (let i = 0; i < planes.length; i++) {
        realTrianglesInstance = clipInstanceAgainstPlane(realTrianglesInstance, planes[i])
        if (realTrianglesInstance == null) {
            return null;
        }
    }
    return realTrianglesInstance;
} 

// instance has vertices and triangles.
function clipInstanceAgainstPlane(instance, plane) {
    let boundingSphere = computeBoundingSphere(instance); 
    let distance = signedDistance(plane, boundingSphere[0]) // distance from object center to plane

    let radius = boundingSphere[1];
    // radius should not be 51.3
    if (distance > radius) {
        return instance; // accepted
    } else if (distance < -radius) {
        return null; // rejected
    } else { // instance intersects the plane.
        let clippedInstance = new Model(instance.vertices, instance.triangles);
        // need to figure out how to handle this data... fix the rendering pipeline...

        // i see, so we call this function 5 times per instance... cant
        // create a helper function to map out the triangles here
        // let triangles = createRealTrianglesList(instance);
        clippedInstance.triangles = clipTrianglesAgainstPlane(clippedInstance.triangles, plane);
        return clippedInstance;
    }

}

function clipTrianglesAgainstPlane(triangles, plane) {
    let clippedTriangles = [];
    for (let i = 0; i < triangles.length; i++) {
        let newTriangles = clipTriangle(triangles[i], plane);
        for (let j = 0; j < newTriangles.length; j++) {
            clippedTriangles.push(newTriangles[j]);
        }
    }

    return clippedTriangles;
}

function clipTriangle(triangle, plane) {
    // i think the best thing to do here is make a helper function to split
    // d0,d1,d2 into two lists. positive and negative.
    // this will help out with the conditionals below.
    let d0 = signedDistance(plane, triangle[0]);
    let d1 = signedDistance(plane, triangle[1]);
    let d2 = signedDistance(plane, triangle[2]);

    let pos = [];
    let neg = [];

    if (d0 >= 0) {
        pos.push(triangle[0]);
    } else {
        neg.push(triangle[0]);
    }
    if (d1 >= 0) {
        pos.push(triangle[1]);
    } else {
        neg.push(triangle[1]);
    }
    if (d2 >= 0) {
        pos.push(triangle[2]);
    } else {
        neg.push(triangle[2]);
    }
    if (pos.length == 3) {
        return [triangle];
    } else if (neg.length == 3) {
        return [];
    } else if (pos.length == 1) {
        let a = pos[0];
        let b = neg[0];
        let c = neg[1];

        let bPrime = computeLinePlaneIntersection(plane, a, b);
        let cPrime = computeLinePlaneIntersection(plane, a, c);

        let triangle1 = [a, bPrime, cPrime, triangle[3]];

        return [triangle1];
    } else if (pos.length == 2) {
        let c = neg[0];
        let a = pos[0];
        let b = pos[1];

        let aPrime = computeLinePlaneIntersection(plane, a, c);
        let bPrime = computeLinePlaneIntersection(plane, b, c);
        let triangle1 = [a, b, aPrime, triangle[3]];
        let triangle2 = [aPrime, b, bPrime, triangle[3]];

        return [triangle1, triangle2];
    }
}

function computeLinePlaneIntersection(plane, a, b) {
    let n = plane.normal;
    let t = (-plane.d - dotProduct(n,a)) / dotProduct(n, subtract(b,a));
    let intersection = add(a, applyScale(subtract(b,a), t));

    return intersection;
}

function createRealTrianglesList(instance) {
    let realTriangles = [];

    let vertices = instance.vertices;
    let triangles = instance.triangles;

    for (let i = 0; i < triangles.length; i ++ ) {
        let triangle = [vertices[triangles[i][0]], vertices[triangles[i][1]], vertices[triangles[i][2]], triangles[i][3]]; 
        realTriangles.push(triangle);
    }
    return realTriangles;
}

function signedDistance(plane, vertex) {
    let normal = plane.normal;

    return (vertex[0] * normal[0]) +
        (vertex[1] * normal[1]) + 
        (vertex[2] * normal[2]) + 
        plane.d;

}

function computeBoundingSphere(instance) {
    // loop through each vertex in the instance.
    // compute avg vertex. this will be the center.
    // loop through again and compute distance from each vertex to center
    // radius is furthest distance..
    //

    let xAvg = 0;
    let yAvg = 0;
    let zAvg = 0;
    let count = 0;
    let vertices = instance.vertices;

    for (count = 0; count < vertices.length; count++) {
        xAvg += vertices[count][0];
        yAvg += vertices[count][1];
        zAvg += vertices[count][2];
    }

    xAvg = xAvg / count;
    yAvg = yAvg / count;
    zAvg = zAvg / count;
    let center = [xAvg, yAvg , zAvg];
    
    let distance = 0;

    for (let i = 0; i < vertices.length; i++) {
        let tempDistance = Math.sqrt(
            Math.pow(xAvg - vertices[i][0], 2) +
            Math.pow(yAvg - vertices[i][1], 2) +
            Math.pow(zAvg - vertices[i][2], 2)
        );
        if (tempDistance > distance) {
            distance = tempDistance;
        }
    }


    return [center, distance];

}

function renderScene() {
    initializeDepthBuffer();
    for (let i = 0; i < scene.length; i++){
        renderInstance(scene[i]);
    }
}

function renderInstance(instance) {
    let projected = [];
    let transformed = [];
    let vertices = instance.model.vertices;
    let triangles = instance.model.triangles;
    for (let i = 0; i < vertices.length; i++) {
        let transformedVertex = applyTransform(vertices[i], instance.transform);
        transformed.push(transformedVertex);
        // so between transformed and projected do clipping.
        // projected.push(projectVertex(transformedVertex));
    }
    let transformedModel = new Model(transformed, triangles);

    let clippedModel = clipInstance(transformedModel, clippingPlanes);
    
    if (clippedModel === null) {
        return;
    }

    //for (let i = 0; i < transformed.length; i++) {
    for (let i = 0; i < clippedModel.vertices.length; i++) {
        projected.push(projectVertex(clippedModel.vertices[i]));
    }

    // doesn't work because it hasn't been projected yet...
    for (let i = 0; i < clippedModel.triangles.length; i++) {
        // so I need to project these vertices
        renderTriangle(clippedModel.triangles[i]);
    }
}

// can merge this into a single transform matrix...
// done by multiplying scale, rotation and translation matrices...
// i don't even feel like this saves time
    // only if I combine them manually without multiplying them each...
function applyTransform(vertex, transform) {
    let scaled = applyScale(vertex, transform.scale);
    let rotated = applyRotation(scaled, transform.rotation);
    let translated = add(rotated, transform.position);
    // return translated;
    return applyCameraTransform(translated);
}

function applyCameraTransform(vertex){
    let inverseCamRotation = applyScale(cameraRotation, -1);
    let inverseCamPosition = applyScale(camera, -1);
    
    let translated = add(vertex, inverseCamPosition);
    let rotated = applyRotation(translated, inverseCamRotation);

    return rotated;
}

function renderTriangle(triangle) {
    // change to drawFilled triangle.. 
    // add depth buffer logic to drawFilledTriangle.


    // so project vertex returns a 2d point. 
    // we SHOULD return a 3d point with the 3rd value containing the triangles un-altered z value
    drawFilledTriangle(projectVertex(triangle[0]),
                          projectVertex(triangle[1]),
                          projectVertex(triangle[2]),
                          triangle[3]);
}

// need to add z point here... that way it can be used by interpolation in drawFilledTriangle
function viewPortToCanvas(x, y, z) {
    return new Point(x * width / v_width, y * height / v_height, z);
}

// make viewProtToCanvas return a 3d point. just keep original z value.
function projectVertex(v) {
    return viewPortToCanvas(v[0] * projection_plane_d / v[2], v[1] * projection_plane_d / v[2], v[2]);

}

function drawWireFrameTriangle(p0, p1, p2, color) {
    drawLine(p0, p1, color);
    drawLine(p1, p2, color);
    drawLine(p2, p0, color);
}

function drawFilledTriangle (p0, p1, p2, color) {
    // so now all points should have an x, y and z value.
    // sort the points so y0 <= y1 <= y2
    if (p1.y < p0.y) { let pT = p1; p1 = p0; p0 = pT; }
    if (p2.y < p0.y) { let pT = p2; p2 = p0; p0 = pT; }
    if (p2.y < p1.y) { let pT = p2; p2 = p1; p1 = pT; }

    // so now all h01 are interpolated using the z value (depth)
    // this means h01 is a mapping of z values corresponding to the points along the x01 line.
    // compute the x coordinates and h values of the triangle edges.
    let x01 = interpolate(p0.y, p0.x, p1.y, p1.x);//xvalues for p0-p1 // ditch h, get z for these values
    let h01 = interpolate(p0.y, p0.z, p1.y, p1.z);//hvalues for p0-p1

    let x12 = interpolate(p1.y, p1.x, p2.y, p2.x);//xvalues for p1-p2
    let h12 = interpolate(p1.y, p1.z, p2.y, p2.z);//h-values for p1-p2

    let x02 = interpolate(p0.y, p0.x, p2.y, p2.x);//xvalues for p0-p2 (long edge)
    let h02 = interpolate(p0.y, p0.z, p2.y, p2.z);//hvalues for p0-p2
    
    // concatenate the short sides.
    x01.pop();
    let x012 = x01.concat(x12);

    h01.pop();
    let h012 = h01.concat(h12); // compute z not h

    // determine which is left and which is right
    let m = Math.floor(x012.length / 2);

    if (x02[m] < x012[m]) {
        var x_left = x02;
        var h_left = h02;

        var x_right = x012;
        var h_right = h012;
    } else {
        var x_left = x012;
        var h_left = h012;

        var x_right = x02;
        var h_right = h02
    }
    
    // compute h value and draw pixel on line between xl and xr.
    for (let y = p0.y; y <= p2.y; y++) {
        let xl = x_left[y - p0.y | 0];
        let xr = x_right[y - p0.y | 0];
        let z_segment = interpolate(xl, h_left[y - p0.y | 0], xr, h_right[y - p0.y | 0]);
        
    
        for (var x = xl; x <= xr; x++) {
            // this line originally should reference the h_segment. see previous chapter code
            let z = z_segment[x-xl | 0];
            let xCanvas = width/2 + (x | 0);
            
            let yCanvas = height/2 - (y | 0) - 1;

            // there was a bug here where it was acessing out of bounds indeces to the debthBuffer.
            // would be better to handle this in my calculations or maybe its a clipping issue so time is not spent checking conditionals here.
            if (xCanvas >= 0 && yCanvas >= 0 && xCanvas < depthBuffer.length && yCanvas < depthBuffer[xCanvas].length && z < depthBuffer[xCanvas][yCanvas]) {

                // now we must assign colors at the triangle level to actually test that this is working.
                putPixel(x, y, color);
                // now to define and initialize depth buffer... 
                // for now we can create a new depth buffer per frame. Maybe there is a more efficient way later on
                // depth buffer needs to account for all triangles.... so it should be initialized at the renderScene level
                // can be a global variable
                // then when renderScene is first called it will be reset. 
                // all values should be set to INFINITY
                depthBuffer[xCanvas][yCanvas] = z;
            }
            
        }
      }
}

// draw a line from point PO to point P1 with color
function drawLine(p0, p1, color) {
    // compute slope.
    let dx = p1.x - p0.x;
    let dy = p1.y - p0.y;
    if (Math.abs(dx) > Math.abs(dy)) {
        //line is horizontal-ish.
        if (dx < 0) {
            // incase p1.x > p0.x
            let pT = p0;
            p0 = p1;
            p1 = pT
        }
        // interpolate with x,y for horizontal function
        let ys = interpolate(p0.x, p0.y, p1.x, p1.y);
        for (let x = p0.x; x <= p1.x; x++) {
            putPixel(x, ys[(x - p0.x | 0)], color);
        }

    } else {
        // Line is vertical-ish
        if (dy < 0) {
            let pT = p0;
            p0 = p1;
            p1 = pT
        }
        // interpolate y,x since line is vertical
        let xs = interpolate(p0.y, p0.x, p1.y, p1.x);
        for (let y = p0.y; y <= p1.y; y++) {
            putPixel(xs[(y - p0.y) | 0], y, color);   
        }
    }
    
    
}

// d = f(i)
// returns array of points along the line
// d is dependent, i is independent
function interpolate(i0, d0, i1, d1) {
    // if i0 == i1, then we just have a point
    if (i0 == i1) {
        return [d0];
      }
    
      var values = [];
      // compute slope
      var a = (d1 - d0) / (i1 - i0);
      var d = d0;
      // add point to values list for each i
      for (var i = i0; i <= i1; i++) {
        values.push(d);
        d += a;
      }

      return values;
    }

// the function that writes a pixes to the canvas buffer
function putPixel(x, y, color) {
    x = width/2 + (x | 0);
    y = height/2 - (y | 0) - 1;
  
    if (x < 0 || x >= width || y < 0 || y >= height) {
        return;
    }
  
    // not sure what canvas pitch is
    let offset = 4*x + canvas_pitch*y;
    canvas_buffer.data[offset++] = color[0];
    canvas_buffer.data[offset++] = color[1];
    canvas_buffer.data[offset++] = color[2];
    canvas_buffer.data[offset++] = 255; // Alpha = 255 (full opacity)
}

// write the canvas buffer to the canvase
function updateCanvas() {
    ctx.putImageData(canvas_buffer, 0, 0);
    resetBuffer();
}

function resetBuffer() {
    for(let i = 0; i < canvas_buffer.data.length; i++){
        canvas_buffer.data[i] = 0;
    }

}

function Point(x, y, z) {
    this.x = x;
    this.y = y;
    this.z = z;
}

function add(v1, v2) {
    return [v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]];
}

function subtract(v1, v2) {
    return [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]];
}

function dotProduct(v1, v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

function applyScale(v1, scalar) {
    return [v1[0] * scalar, v1[1] * scalar, v1[2] * scalar]
}

function applyRotation(v1, rotations) {
    let rotationMatrix = constructRotationMatrix(rotations);
    return matMult(rotationMatrix, v1);
}

function invert(v1){
    let inverted = [0, 0, 0];
    if (v1[0] != 0) {
        inverted[0] = 1/v1[0];
    }
    if (v1[1] != 0) {
        inverted[1] = 1/v1[1];
    }
    if (v1[2] != 0) {
        inverted[2] = 1/v1[2];
    }
    return inverted;
}

function constructRotationMatrix(rotations) {
    let roll = rotations[0] * Math.PI / 180.0;  // x
    let pitch = rotations[1] * Math.PI / 180.0;   // y
    let yaw = rotations[2] * Math.PI / 180.0;    // z

    let sinPitch = Math.sin(pitch);
    let cosPitch = Math.cos(pitch);
    
    let sinRoll = Math.sin(roll);
    let cosRoll = Math.cos(roll);

    let sinYaw = Math.sin(yaw);
    let cosYaw = Math.cos(yaw);

    let rotationMatrix = [
        [cosYaw*cosPitch, cosYaw*sinPitch*sinRoll - sinYaw*cosRoll, cosYaw*sinPitch*cosRoll + sinYaw*sinRoll],
        [sinYaw*cosPitch, sinYaw*sinPitch*sinRoll + cosYaw*cosRoll, sinYaw*sinPitch*cosRoll - cosYaw*sinRoll],
        [-sinPitch, cosPitch*sinRoll, cosPitch*cosRoll]
    ];

    return rotationMatrix;



}

// this is a junky implementation of matrix multiplication, will need to fix this. better for 4x4 anyways to account for translations and scales....
function matMult(m1, m2) {
    // assuming 3x3, 3x1
    let m3 = [
        dotProduct(m1[0], m2),
        dotProduct(m1[1], m2),
        dotProduct(m1[2], m2),
    ];
    return m3;
}

function computeCamMovement(angle) {
    let theta = angle * Math.PI / 180.0;
    let zMovement = Math.cos(theta);
    let xMovement = Math.sin(theta);
    return [zMovement, xMovement];
}