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

let black = [0,0,0];
let green = [0,255,0];
let blue = [0,0,255];
let red = [255,0,0];


// The four "front vertices"
let vAf = [-2, -0.5, 5];
let vBf = [-2,  0.5, 5];
let vCf = [-1,  0.5, 5];
let vDf = [-1, -0.5, 5];

// The four "back" vertices
let vAb = [-2, -0.5, 6];
let vBb = [-2,  0.5, 6];
let vCb = [-1,  0.5, 6];
let vDb = [-1, -0.5, 6];

// front face
drawLine(projectVertex(vAf), projectVertex(vBf), blue);
drawLine(projectVertex(vBf), projectVertex(vCf), blue);
drawLine(projectVertex(vCf), projectVertex(vDf), blue);
drawLine(projectVertex(vDf), projectVertex(vAf), blue);

// The back face
drawLine(projectVertex(vAb), projectVertex(vBb), red);
drawLine(projectVertex(vBb), projectVertex(vCb), red);
drawLine(projectVertex(vCb), projectVertex(vDb), red);
drawLine(projectVertex(vDb), projectVertex(vAb), red);

// The front-to-back edges
drawLine(projectVertex(vAf), projectVertex(vAb), green);
drawLine(projectVertex(vBf), projectVertex(vBb), green);
drawLine(projectVertex(vCf), projectVertex(vCb), green);
drawLine(projectVertex(vDf), projectVertex(vDb), green);


updateCanvas();


function viewPortToCanvas(x, y) {
    return new Point(x * width / v_width, y * height / v_height);
}

function projectVertex(v) {
    return viewPortToCanvas(v[0] * projection_plane_d / v[2], v[1] * projection_plane_d / v[2]);

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
}

function Point(x, y) {
    this.x = x;
    this.y = y;
}
