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
    [4, 0, 3, red],
    [4, 3, 7, red],
    [5, 4, 7, red],
    [5, 7, 6, red],
    [1, 5, 6, red],
    [1, 6, 2, red],
    [4, 5, 1, red],
    [4, 1, 0, red],
    [2, 6, 7, red],
    [2, 7, 3, red]
]


let cubeModel = new Model(cubeVertices, cubeTriangles);

let cube1 = new Instance(cubeModel, new Transform(1, [0, 0, 0], [-1.5, 0, 7]));
let cube3 = new Instance(cubeModel, new Transform(1, [0, 0, 0], [-1.5, 0, -7]));
let cube2 = new Instance(cubeModel, new Transform(1, [45, 45, 45], [1.25, 2, 10.5]));

let scene = [cube1];
let s2 = Math.sqrt(2);
let planeNear = new Plane([0, 0, 1], -projection_plane_d);
let planeLeft = new Plane([s2, 0, s2], 0);
let planeRight = new Plane([-s2, 0, s2], 0);
let planeTop = new Plane([0, -s2, s2], 0);
let planeBottom = new Plane([0, s2, s2], 0);

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

function Plane(normal, d) {
    this.normal = normal;
    this.d = d;
}

// function clipScene(scene, planes) {
//     let clippedInstances = [];
//     for (let i = 0; i < scene.length; i++){
//         // to-do: define planes (compute them) and implement clipInstance
//         let clipped = clipInstance(scene[i], planes);
//         if (clipped != null){
//             clippedInstances.push(clipped);
//         }
//     }
//     return clippedInstances;
// }

// need to define planes.
// planes are defined as.....
// A, B, C, D

function clipInstance(instance, planes) {
    console.log("START CLIPPING NEW INSTANCE");
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
        console.log("accept whole object");
        return instance; // accepted
    // } else if (distance < -radius) {
    //     console.log("reject whole object"); // i guess this is somehow a problem.
    //     return null; // rejected
    } else { // instance intersects the plane.
        let clippedInstance = new Model(instance.vertices, instance.triangles);
        console.log("clip triangles");
        // need to figure out how to handle this data... fix the rendering pipeline...

        // i see, so we call this function 5 times per instance... cant
        // create a helper function to map out the triangles here
        // let triangles = createRealTrianglesList(instance);
        clippedInstance.triangles = clipTrianglesAgainstPlane(clippedInstance.triangles, plane);
        console.log(clippedInstance.triangles.length);
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
    console.log("clipping new triange");
    // i think the best thing to do here is make a helper function to split
    // d0,d1,d2 into two lists. positive and negative.
    // this will help out with the conditionals below.
    let d0 = signedDistance(plane, triangle[0]);
    let d1 = signedDistance(plane, triangle[1]);
    let d2 = signedDistance(plane, triangle[2]);

    let pos = [];
    let neg = [];

    if (d0 >= 0) {
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
;;;;;;;;;;;;;;;;;;;;;;;;
    if (pos.length == 3) {
        console.log("keep original triangle");
        return [triangle];
    } else if (neg.length == 3) {
        console.log("reject entire triangle");
        return []; // is this the right thing to do here? I don't think so.
        // maybe better to return null and then add a null check in clipTrainglesAgainstPlane.
        // this is actually correct because this function will return a list with 0 to 2 triangles.
    } else if (pos.length == 1) {
        console.log("clip triangle into one triangle");
        // how to handle this data so that i can get the point coresponding with the distance.... oh just push point and not the distance.

        // another problem... as of now this triangle will contain a mapping to the vertices, not the actual vertices themselves.....
        // i.e a triangle might be [0,1,2] signifying vertexes 0, 1 and 2 instead of holding the actual value of these vertices.
        let a = pos[0];
        let b = neg[0];
        let c = neg[1];

        // add plane intersection function here.
        let bPrime = computeLinePlaneIntersection(plane, a, b);
        let cPrime = computeLinePlaneIntersection(plane, a, c);

        let triangle1 = [a, bPrime, cPrime, triangle[3]];

        return [triangle1];
    } else if (pos.length == 2) {
        console.log("clip triangle into two triangles");
        let c = neg[0];
        let a = pos[0];
        let b = pos[1];

        let aPrime = computeLinePlaneIntersection(plane, a, c);
        let bPrime = computeLinePlaneIntersection(plane, b, c);
        let triangle1 = [a, b, aPrime, triangle[3]];// triangle 3 color
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
        let triangle = [vertices[triangles[i][0]], vertices[triangles[i][1]], vertices[triangles[i][2]], triangles[i][3]]; // will this have a problem with previous clipping? i dont think so.
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

    // problem with this is that this should occur after model translations and camera translations have happened and rotations
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
    yavg = yAvg / count;
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
    // add clipping here. // we're going to skip "clipScene" function and go straight to clipInstance
    let transformedModel = new Model(transformed, triangles);

    let clippedModel = clipInstance(transformedModel, clippingPlanes);
    
    if (clippedModel === null) {
        return;
    }

    //for (let i = 0; i < transformed.length; i++) {
    for (let i = 0; i < clippedModel.vertices.length; i++) {
        projected.push(projectVertex(clippedModel.vertices[i]));
        // projected.push(projectVertex(transformed[i]));
    }

    // casnt render clippedModel.triangles like this because I've changed the data structure of the triangles
    // clippedModel has real triangles, while renderTriangleFunction is expecting a index mapping to vertices for the triangles.
    // for (let i = 0; i < triangles.length; i++) {
    //     console.log("num triangles rendering: ", triangles.length);
    //     // so in theory I could just pass in triangles here and modify the function.
    //     renderTriangle(triangles[i], projected);
    // }

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

function renderTriangle(triangle, projected) {
    drawWireFrameTriangle(projected[triangle[0]],
                          projected[triangle[1]],
                          projected[triangle[2]],
                          triangle[3]);
}

function renderTriangle(triangle) {
    drawWireFrameTriangle(projectVertex(triangle[0]),
                          projectVertex(triangle[1]),
                          projectVertex(triangle[2]),
                          triangle[3]);
}

function viewPortToCanvas(x, y) {
    return new Point(x * width / v_width, y * height / v_height);
}

function projectVertex(v) {
    return viewPortToCanvas(v[0] * projection_plane_d / v[2], v[1] * projection_plane_d / v[2]);

}

function drawWireFrameTriangle(p0, p1, p2, color) {
    drawLine(p0, p1, color);
    drawLine(p1, p2, color);
    drawLine(p2, p0, color);
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

function Point(x, y) {
    this.x = x;
    this.y = y;
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
    let theta = angle * Math.PI / 180.0;;
    let zMovement = Math.cos(theta);
    let xMovement = Math.sin(theta);
    return [zMovement, xMovement];
}