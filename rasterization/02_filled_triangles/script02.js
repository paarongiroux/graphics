let c = document.getElementById("myCanvas");
let ctx = c.getContext("2d");
let width = c.width;
let height = c.height;
let canvas_buffer = ctx.getImageData(0, 0, width, height);
let canvas_pitch = canvas_buffer.width * 4;

let black = [0,0,0];
let green = [255,255,0];
let tri0 = new Point(-200, -250, 0.5);
let tri1 = new Point(200, -50, 0);
let tri2 = new Point(20, 250, 1);

drawFilledTriangle(tri0, tri1, tri2, green);
drawWireFrameTriangle(tri0, tri1, tri2, black);

updateCanvas();

function drawFilledTriangle (p0, p1, p2, color) {
    // sort the points so y0 <= y1 <= y2
    if (p1.y < p0.y) { let pT = p1; p1 = p0; p0 = pT; }
    if (p2.y < p0.y) { let pT = p2; p2 = p0; p0 = pT; }
    if (p2.y < p1.y) { let pT = p2; p2 = p1; p1 = pT; }

    // compute the x coordinates and h values of the triangle edges.
    let x01 = interpolate(p0.y, p0.x, p1.y, p1.x);
    let h01 = interpolate(p0.y, p0.h, p1.y, p1.h);

    let x12 = interpolate(p1.y, p1.x, p2.y, p2.x);
    let h12 = interpolate(p1.y, p1.h, p2.y, p2.h);

    let x02 = interpolate(p0.y, p0.x, p2.y, p2.x); //long edge
    let h02 = interpolate(p0.y, p0.h, p2.y, p2.h);
    
    // concatenate the short sides.
    x01.pop();
    let x012 = x01.concat(x12);

    h01.pop();
    let h012 = h01.concat(h12);

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
        let xl = x_left[y - p0.y] | 0;
        let xr = x_right[y - p0.y] | 0;
        let h_segment = interpolate(xl, h_left[y - p0.y], xr, h_right[y - p0.y]);
    
        for (var x = xl; x <= xr; x++) {
          putPixel(x, y, multiply(h_segment[x - xl], color));
        }
      }
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
}

function dot_product(p1, p2) {
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

function multiply(k, vec) {
    return [k*vec[0], k*vec[1], k*vec[2]];
  }

function Point(x, y, h) {
    this.x = x;
    this.y = y;
    this.h = h;
}