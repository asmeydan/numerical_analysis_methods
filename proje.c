#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

// math serbest mi
// polinom haricindeki fonksiyonlar nasıl olacak

typedef struct
{
    int type;
    double coefficient;
    double degree;
} POLY;

double parse_number(char *expr, int *index);
double parse_factor(char *expr, int *index, int x);
double parse_term(char *expr, int *index, int x);
double parse_expression(char *expr, int *index, int x);
double parse_exponent(char *expr, int *index, int x);
void newtonRaphsonMethod();
void regulaFalsiMethod();
void bisectionMethod();
void printEquation(POLY *equation, int equationCount);
POLY *getEquation(int equationCount);
POLY *getDerivative(POLY *equation, int equationCount);
double calcEquation(POLY *equation, int equationCount, double x);
int menu();

int main()
{

    char expression[256];
    printf("Eşitsizlik gir: ");
    scanf("%s", expression);

    double x = 2.0; // Example value for x
    int index = 0;
    double result = parse_expression(expression, &index, x);
    printf("Result of '%s' with x=%.2lf: %lf\n", expression, x, result);

    int choice;

    choice = menu();

    switch (choice)
    {
    case 1:
        bisectionMethod();
        break;

    case 2:
        regulaFalsiMethod();
        break;

    case 3:
        newtonRaphsonMethod();
        break;

    default:
        break;
    }

    return 0;
}

void newtonRaphsonMethod()
{
    double start, x, error, expectedError, realValue, result, tmp, tmp2;
    int i, equationCount, iterMax, iter = 0, errorType;

    printf("\nenter element count: ");
    scanf("%d", &equationCount);

    POLY *equation = getEquation(equationCount);
    printEquation(equation, equationCount);
    POLY *derivative = getDerivative(equation, equationCount);
    printEquation(derivative, equationCount);

    printf("Enter a max iteration value: ");
    scanf("%d", &iterMax);
    printf("Enter x value: ");
    scanf("%lf", &x);

    do
    {
        printf("0- real value\n1- otherwise \nerror type: ");
        scanf("%d", &errorType);
    } while (errorType < 0 || errorType > 1);
    if (errorType == 0)
    {
        printf("real value: ");
        scanf("%lf", &realValue);
    }
    printf("expected error: ");
    scanf("%lf", &expectedError);

    do
    {
        iter++;
        // PROBLEM
        tmp = calcEquation(equation, equationCount, x) / calcEquation(derivative, equationCount, x);
        tmp2 = x;
        x = x - tmp;
        result = calcEquation(equation, equationCount, x);
        printf("tmp: %lf  mid: %lf  result: %lf\n", tmp, x, result);

        if (errorType == 0)
        {
            error = fabs(realValue - result);
        }
        else if (errorType == 1)
        {
            error = fabs(x - tmp2);
        }
        printf("Regula Falsi iter: %d  x: %lf  result: %lf  error: %lf\n", iter, x, result, error);
    } while (error > expectedError && iter < iterMax);
    printf("result: %lf\n", x);

    free(equation);
    free(derivative);
}

void regulaFalsiMethod()
{
    double start, left, right, mid, error, expectedError, realValue, result, tmp;
    int i, equationCount, iterMax, iter = 0, errorType;

    printf("\nenter element count: ");
    scanf("%d", &equationCount);

    POLY *equation = getEquation(equationCount);
    printEquation(equation, equationCount);

    printf("Enter a max iteration value: ");
    scanf("%d", &iterMax);
    printf("Enter left value for [left,right]: ");
    scanf("%lf", &left);
    printf("Enter right value for [left,right]: ");
    scanf("%lf", &right);

    do
    {
        printf("0- real value\n1- otherwise \nerror type: ");
        scanf("%d", &errorType);
    } while (errorType < 0 || errorType > 1);
    if (errorType == 0)
    {
        printf("real value: ");
        scanf("%lf", &realValue);
    }
    printf("expected error: ");
    scanf("%lf", &expectedError);

    if (calcEquation(equation, equationCount, left) * calcEquation(equation, equationCount, right) > 0)
    {
        printf("invalid range ");
    }
    else
    {
        do
        {
            iter++;
            // PROBLEM
            tmp = calcEquation(equation, equationCount, right) - calcEquation(equation, equationCount, left);
            mid = left - (calcEquation(equation, equationCount, left) * ((right - left) / tmp));
            result = calcEquation(equation, equationCount, mid);
            // printf("tmp: %lf  mid: %lf  result: %lf\n", tmp, mid, result);
            if (calcEquation(equation, equationCount, mid) * calcEquation(equation, equationCount, left) < 0)
            {
                right = mid;
            }
            else if (calcEquation(equation, equationCount, mid) * calcEquation(equation, equationCount, right) < 0)
            {
                left = mid;
            }

            if (errorType == 0)
            {
                error = fabs(realValue - result);
            }
            else if (errorType == 1)
            {
                error = (right - left) / pow(2, iter);
            }
            printf("Regula Falsi iter: %d  x: %lf  result: %lf  error: %lf\n", iter, mid, result, error);
        } while (error > expectedError && iter < iterMax);
        printf("result: %lf\n", mid);
    }

    free(equation);
}

void bisectionMethod()
{
    double start, left, right, mid, error, expectedError, realValue, result;
    int i, equationCount, iterMax, iter = 0, errorType;

    printf("\nenter element count: ");
    scanf("%d", &equationCount);

    POLY *equation = getEquation(equationCount);
    printEquation(equation, equationCount);

    printf("Enter a max iteration value: ");
    scanf("%d", &iterMax);
    printf("Enter left value for [left,right]: ");
    scanf("%lf", &left);
    printf("Enter right value for [left,right]: ");
    scanf("%lf", &right);

    do
    {
        printf("0- real value\n1- otherwise \nerror type: ");
        scanf("%d", &errorType);
    } while (errorType < 0 || errorType > 1);
    if (errorType == 0)
    {
        printf("real value: ");
        scanf("%lf", &realValue);
    }
    printf("expected error: ");
    scanf("%lf", &expectedError);

    if (calcEquation(equation, equationCount, left) * calcEquation(equation, equationCount, right) > 0)
    {
        printf("invalid range ");
    }
    else
    {
        do
        {
            iter++;
            mid = (left + right) / 2;
            result = calcEquation(equation, equationCount, mid);
            if (calcEquation(equation, equationCount, mid) * calcEquation(equation, equationCount, left) < 0)
            {
                right = mid;
            }
            else if (calcEquation(equation, equationCount, mid) * calcEquation(equation, equationCount, right) < 0)
            {
                left = mid;
            }

            if (errorType == 0)
            {
                error = fabs(realValue - result);
            }
            else if (errorType == 1)
            {
                error = (right - left) / pow(2, iter);
            }
            printf("Bisection iter: %d  x: %lf  result: %lf  error: %lf\n", iter, mid, result, error);
        } while (error > expectedError && iter < iterMax);
        printf("result: %lf\n", mid);
    }

    free(equation);
}

POLY *getDerivative(POLY *equation, int equationCount)
{
    int i;
    double tmp, result = 0;
    POLY *derivative = (POLY *)malloc(equationCount * sizeof(POLY));

    for (i = 0; i < equationCount; i++)
    {
        if (equation[i].type == 1)
        {
            derivative[i].type = 1;
            derivative[i].coefficient = equation[i].coefficient * equation[i].degree;
            derivative[i].degree = equation[i].degree - 1;
        }
    }

    return derivative;
}

double calcEquation(POLY *equation, int equationCount, double x)
{
    int i;
    double tmp, result = 0;

    for (i = 0; i < equationCount; i++)
    {
        if (equation[i].type == 1)
        {
            tmp = pow(x, equation[i].degree);
            tmp *= equation[i].coefficient;
        }
        result += tmp;
    }

    return result;
}

void printEquation(POLY *equation, int equationCount)
{
    int i;

    for (i = 0; i < equationCount - 1; i++)
    {
        if (equation[i].type == 1)
        {
            printf("%lf x^%lf + ", equation[i].coefficient, equation[i].degree);
        }
    }

    if (equation[i].type == 1)
    {
        printf("%lf x^%lf", equation[i].coefficient, equation[i].degree);
    }
    printf("\n");
}

POLY *getEquation(int equationCount)
{
    POLY *res = (POLY *)malloc(equationCount * sizeof(POLY));
    int i, type;

    for (i = 0; i < equationCount; i++)
    {
        printf("\n1- polynomial\n");
        printf("2- exponential\n");
        printf("3- trigonometric\n");
        printf("4- logarithmic\n");
        printf("5- inverse trigonometric\n");
        do
        {
            printf("type: ");
            scanf("%d", &type);
        } while (type <= 0 || type > 5);
        res[i].type = type;
        if (res[i].type == 1)
        {
            printf("\nEnter x's degree: ");
            scanf("%lf", &res[i].degree);
            printf("Enter x's coefficient: ");
            scanf("%lf", &res[i].coefficient);
        }
    }

    return res;
}

double parse_exponent(char *expr, int *index, int x) {
    double base = parse_factor(expr, index, x);
    
    while (expr[*index] == '^') {
        (*index)++;  // Skip '^'
        double exponent = parse_factor(expr, index, x);
        base = pow(base, exponent);
    }
    
    return base;
}

double parse_number(char *expr, int *index)
{
    double number = 0;
    sscanf(expr + *index, "%lf", &number);
    while (isdigit(expr[*index]) || expr[*index] == '.')
    {
        (*index)++;
    }
    return number;
}

double parse_factor(char *expr, int *index, int x)
{
    double result = 0.0;

    if (expr[*index] == '(')
    {
        (*index)++;
        result = parse_expression(expr, index, x);
        (*index)++;
    }
    else if (strncmp(expr + *index, "sin", 3) == 0)
    {
        (*index) += 3;
        (*index)++;
        result = sin(parse_expression(expr, index, x));
        (*index)++;
    }
    else if (strncmp(expr + *index, "log_", 4) == 0)
    {
        (*index) += 4;
        double base = parse_expression(expr, index, x);
        (*index)++;
        result = log(parse_expression(expr, index, x)) / log(base);
        (*index)++;
    }
    else if (strncmp(expr + *index, "e^", 2) == 0)
    {
        (*index) += 2;
        result = exp(parse_expression(expr, index, x));
    }
    else if (strncmp(expr + *index, "x", 1) == 0)
    {
        (*index) += 1;
        result = x;
    }
    else
    {
        result = parse_number(expr, index);
    }

    return result;
}

double parse_term(char *expr, int *index, int x)
{
    double result = parse_exponent(expr, index, x);

    while (expr[*index] == '*' || expr[*index] == '/')
    {
        char op = expr[*index];
        (*index)++;
        double next_factor = parse_exponent(expr, index, x);

        if (op == '*')
        {
            result *= next_factor;
        }
        else
        {
            result /= next_factor;
        }
    }

    return result;
}

double parse_expression(char *expr, int *index, int x)
{
    double result = parse_term(expr, index, x);

    while (expr[*index] == '+' || expr[*index] == '-')
    {
        char op = expr[*index];
        (*index)++;
        double next_term = parse_term(expr, index, x);

        if (op == '+')
        {
            result += next_term;
        }
        else
        {
            result -= next_term;
        }
    }

    return result;
}

int menu()
{
    int method, elementCount;
    printf("\nBisection Method:1\n");
    printf("Regula-Falsi Method:2 \n");
    printf("Newton-Raphson Method:3 \n");
    printf("Inverse of a matrix:4\n");
    printf("Gauss-Jordan elimination:5\n");
    printf("Gauss seidel Method:6\n");
    printf("Numerical Differentiation:7\n");
    printf("Simpson Method:8\n");
    printf("Trapezoidal Method:9\n");
    printf("Gregory Newton Interpolation:10\n");
    printf("Method of choice: ");
    scanf("%d", &method);
    return method;
}