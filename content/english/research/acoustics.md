---
title: "concert hall simulation"
mathjax: true
katex: true
markup: 'mmark'
math: true
draft: false
---

This simulation attempts to recreate the distribution of sound in square concert halls. The set-up is seen in the figure below, where the radii are calculate through geometry

![setup](/acoustics_setup.svg#center?w100)

The radii come out to be 

$$\begin{aligned}
r_1 &= d_1 \sqrt{\left(\frac{y}{2d_1-x}\right)^2 + 1} 
\\\\
r_3 &= d_3 \sqrt{\left(\frac{y}{2d_1+x}\right)^2 + 1} 
\\\\
r_2 &= \sqrt{\left(\frac{y(2d_1-x)}{d_1}\right)^2 + (d_1 - x)^2} 
\\\\
r_4 &= \sqrt{\left(\frac{y(2d_1+x)}{d_1}\right)^2 + (d_1 + x)^2} 
\\\\
r_5 &= \sqrt{x^2 + y^2}
\end{aligned}$$

Then, after finding the sound level through 
$$\beta = \log{\left(\frac{I}{I_0}\right)}$$
where $I$ is the intensity and $A$ is amplitude of each (total) sound wave. An example result of the code is this image!

![BSH](/acoustics_BSH.png#center)

An absorbtion can also be simply applied through the absorption coefficienct $\alpha = I_{\alpha} / I$ which results in this:

![BSH](/acoustics_smooth.png#center)