---
title: "code and graphic design"
# date: 2022-12-20
date: 2022-12-20T15:17:05-08:00
tags: ["python"]
---
I wanted to look into the intersection between graphic design and code.
I found a program called Processing, where you can code in Java, C++ or Python.
After trying a few simple things, this animation seemed to be pretty cool:

# {{< video src="/out.mp4" type="video/mp4" preload="auto">}}

The code for it is below. There are two functions that automatically run, 
`setup()` (runs once like `__init__`) and `draw()` (continuously runs). So I just
continuously increased the radius at which that ellipse is at, as well as it's angle with
a couple sines and cosines to make it dynamic.

```
import math
import random

class DrawEllipse():
    def __init__(self):
        self.angle = 0
        self.r = 1
        
    def draw_circle(self, r=50):
    	# increase radius and angle incrementally
        self.angle += 0.05
        self.r += 0.1
        # get cartesian coords
        x_coord = r * math.cos(self.angle)
        y_coord = r * math.sin(self.angle)
        # color and draw ellipse
        fill(10, 50, 150)
        ellipseMode(CENTER)
        ellipse(x_coord, y_coord, self.r * math.cos(self.angle), self.r * math.sin(self.angle))
        
    def get_angle(self):
        return self.angle
        
d = DrawEllipse()

def setup():
    size(300, 300)
    background(20, 20, 20)
    
def draw():
    translate(width / 2, height / 2)
    d.draw_circle()
```
