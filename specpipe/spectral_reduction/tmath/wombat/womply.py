def womply(a,i,x):
    """poly evaluator for womashrebin
;
;             I+1
;  EVALUATES  SUM A(J) * X**(I+1-J)
;             J=1
;
;  I IS THUS THE ORDER OF THE POLY, WHICH HAS I+1 COEFFICIENTS,
;  AND A(I+1) IS THE CONSTANT TERM.
;

    x is not array-like

    """
    if (i < 0):
        polyval=0.
        print('Value of I out of bounds in womply')
        return polyval
    polyval=a[0]
    if (i == 0):
        return polyval
    for j in range(2,i+2):
        polyval=x*polyval+a[j-1]
    return polyval

