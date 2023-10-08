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
        cameraRotation[1] -= 2;
        renderScene();
        updateCanvas();
    }
    if (code==="KeyL"){
        cameraRotation[1] += 2;
        console.log(cameraRotation[1]);
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
    
    console.log(`Key pressed ${name} \r\n Key code value: ${code}`);
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
let cube2 = new Instance(cubeModel, new Transform(1, [45, 45, 45], [1.25, 2, 10.5]));

let scene = [cube1, cube2];



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

function renderScene() {
    for (let i = 0; i < scene.length; i++){
        renderInstance(scene[i]);
    }
}

function renderScene2() {
    // to-do: create this function
    // to-do: define camera object with position and orientation
    mCamera = createCameraMatrix(camera.position, camera.orientation);
    for (let i = 0; i < scene.length; i++){
        // create a 4x4matmult function
        // define transform field for instance type. i guess this will be its own matrix??
        let mTrans = matMult(mCamera, scene[i].transform);
        // create renderModel function
        renderModel(scene[i].model, mTrans);
    }
}

function renderModel(model, transform) {
    let projected = [];
    for (let i = 0; i < vertices.length; i++) {
        let transformedVertex = matmult(transform, vertices[i])
        projected.push(projectVertex(transformedVertex));
    }
    for (let i = 0; i < triangles.length; i++) {
        renderTriangle(triangles[i], projected);
    }
}

function renderInstance(instance) {
    let projected = [];
    let vertices = instance.model.vertices;
    let triangles = instance.model.triangles;
    for (let i = 0; i < vertices.length; i++) {
        let transformedVertex = applyTransform(vertices[i], instance.transform);
        projected.push(projectVertex(transformedVertex));
    }
    for (let i = 0; i < triangles.length; i++) {
        renderTriangle(triangles[i], projected);
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