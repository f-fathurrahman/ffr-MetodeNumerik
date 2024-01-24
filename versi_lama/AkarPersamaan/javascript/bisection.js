function bisection(f, x1, x2)
{
    f1 = f(x1);
    f2 = f(x2);

    console.log("f1 = ", f1);
    console.log("f2 = ", f2);
    if(f1*f2 > 0.0) {
        console.log("Tanda f1 dan f2 sama");
        return;
    }

    xr = 0.5*(x1 + x2);
    console.log("xr = ", xr);


}

function my_func1(x)
{
    return Math.pow(x, 3) - 35;
}

function test_bisection()
{
    bisection(my_func1, 1.0, 4.0);
}

test_bisection()