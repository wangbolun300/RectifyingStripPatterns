This is for you to test the initialization -> optimization -> binormal strip optimization -> rectifying developable optimization.

"Wave.obj" is the triangle mesh,
"saved_traced_*" are the traced curves of \theta = 60. You can load them one by one. 
"smooth.csv" is a random smooth scalar field without singularities. Load it as the start of initialization.
"initial_ls.csv" is the scalar field generally follows the traced P-curves.
"60.csv" and "60_b.csv" are the scalar field and the binormal vector field of \theta = 60.
"ply*" are the strips extracted from the fields.
"opt*" are the optimized high-quality binormal strips optimized from "ply*".
"crease*" are the rectifying developables optimized from "opt*".
