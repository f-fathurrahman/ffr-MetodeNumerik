To use font TeX font instead of STIX:

/home/efefer/miniconda3/lib/python3.7/site-packages/notebook

MathJax.Hub.Config({
        "HTML-CSS": {
            availableFonts: ["TeX"],
            preferredFonts: "TeX",
         }
})


Unzip the MathJax file and go into the folder

copy jax/output/HTML-CSS/fonts/TeX into directoy ../notebook/static/components/MathJax/jax/output/HTML-CSS/fonts/

copy fonts/HTML-CSS/TeX into ../notebook/static/components/MathJax/fonts/HTML-CSS/

open ../notebook/static/notebook/js/main.min.js, search for availableFonts. It should be around line 14894. Change it to

...
availableFonts: ["STIX-Web","TeX"],
imageFont: null;
preferredFont: "TeX",
webFont: "TeX",
...


Finite difference coefficients are taken from:
http://web.media.mit.edu/~crtaylor/calculator.html
