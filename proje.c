#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

typedef struct
{
    int type;
    double coefficient;
    double degree;
} POLY;

void gregoryNewtonMethod();
int factorial(int n);
void gaussSeidel();
void gaussElimination();
void printMatrix(double **matrix, int rows, int cols);
void inverseMatrix();
void trapezoidalMethod();
void simpsonMethod();
double differentiation(char *expression, int differentiationType, double h, double x);
void numericDerivative();
double parse_number(char *expr, int *index);
double parse_factor(char *expr, int *index, double x);
double parse_term(char *expr, int *index, double x);
double parse_expression(char *expr, int *index, double x);
double parse_exponent(char *expr, int *index, double x);
void printExpression(char *expression);
double calcExpression(char *expression, double x);
char *getExpression();
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

    case 4:
        inverseMatrix();
        break;

    case 5:
        gaussElimination();
        break;

    case 6:
        gaussSeidel();
        break;

    case 7:
        numericDerivative();
        break;

    case 8:
        simpsonMethod();
        break;

    case 9:
        trapezoidalMethod();
        break;

    case 10:
        gregoryNewtonMethod();
        break;

    default:
        break;
    }

    return 0;
}

void gregoryNewtonMethod()
{
    int n, i, j;
    double x, result = 0;

    printf("Veri nokta sayısını girin (n): ");
    scanf("%d", &n);

    double *X = (double *)malloc(n * sizeof(double));
    double *Y = (double *)malloc(n * sizeof(double));
    double **differenceTable = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) {
        differenceTable[i] = (double *)malloc(n * sizeof(double));
    }

    printf("X ve Y değerlerini girin:\n");
    for (i = 0; i < n; i++) {
        printf("X[%d]: ", i);
        scanf("%lf", &X[i]);
        printf("Y[%d]: ", i);
        scanf("%lf", &Y[i]);
        differenceTable[i][0] = Y[i];
    }

    printf("Enterpolasyon yapılacak x değerini girin: ");
    scanf("%lf", &x);

    for (j = 1; j < n; j++) {
        for (i = 0; i < n - j; i++) {
            differenceTable[i][j] = differenceTable[i + 1][j - 1] - differenceTable[i][j - 1];
        }
    }

    double term = 1;
    result = Y[0];
    for (i = 1; i < n; i++) {
        term *= (x - X[i - 1]);
        result += (term * differenceTable[0][i]) / (double)factorial(i);
    }

    printf("Enterpolasyon sonucu f(%lf) = %lf\n", x, result);

    for (i = 0; i < n; i++) {
        free(differenceTable[i]);
    }
    free(differenceTable);
    free(X);
    free(Y);
}

int factorial(int n)
{
    if (n == 0 || n == 1)
        return 1;
    else
        return n * factorial(n - 1);
}

void gaussSeidel() {
    int N, i, j, k, iterMax;
    double tol;

    printf("Matrisin boyutunu girin (N): ");
    scanf("%d", &N);

    double **matrix = (double **)malloc(N * sizeof(double *));
    double *results = (double *)malloc(N * sizeof(double));
    double *solutions = (double *)malloc(N * sizeof(double));
    double *oldSolutions = (double *)malloc(N * sizeof(double));

    for (i = 0; i < N; i++) {
        matrix[i] = (double *)malloc((N + 1) * sizeof(double));
    }

    printf("Matrisin katsayılarını girin:\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("Matris[%d][%d]: ", i, j);
            scanf("%lf", &matrix[i][j]);
        }
    }

    printf("Sonuçlar vektörünü girin:\n");
    for (i = 0; i < N; i++) {
        printf("Sonuç[%d]: ", i);
        scanf("%lf", &results[i]);
        matrix[i][N] = results[i];
    }

    printf("Maksimum iterasyon sayısını girin: ");
    scanf("%d", &iterMax);
    printf("Hata toleransını girin: ");
    scanf("%lf", &tol);

    for (i = 0; i < N; i++) {
        solutions[i] = 0.0; // İlk tahmin olarak sıfırdan başlıyoruz
    }

    for (k = 0; k < iterMax; k++) {
        for (i = 0; i < N; i++) {
            oldSolutions[i] = solutions[i];
        }

        for (i = 0; i < N; i++) {
            double sum = results[i];
            for (j = 0; j < N; j++) {
                if (j != i) {
                    sum -= matrix[i][j] * solutions[j];
                }
            }
            solutions[i] = sum / matrix[i][i];
        }

        // Hata kontrolü
        double maxError = 0.0;
        for (i = 0; i < N; i++) {
            double error = fabs(solutions[i] - oldSolutions[i]);
            if (error > maxError) {
                maxError = error;
            }
        }

        if (maxError < tol) {
            break;
        }
    }

    printf("Çözümler:\n");
    for (i = 0; i < N; i++) {
        printf("x%d = %lf\n", i + 1, solutions[i]);
    }

    for (i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(results);
    free(solutions);
    free(oldSolutions);
}

void gaussElimination() {
    int N;
    int i, j, k;

    printf("Matrisin boyutunu girin (N): ");
    scanf("%d", &N);

    double **matrix = (double **)malloc(N * sizeof(double *));
    double *results = (double *)malloc(N * sizeof(double));
    double *solutions = (double *)malloc(N * sizeof(double));
    for (i = 0; i < N; i++) {
        matrix[i] = (double *)malloc((N + 1) * sizeof(double));
    }

    printf("Matrisin katsayılarını girin:\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("Matris[%d][%d]: ", i, j);
            scanf("%lf", &matrix[i][j]);
        }
    }

    printf("Sonuçlar vektörünü girin:\n");
    for (i = 0; i < N; i++) {
        printf("Sonuç[%d]: ", i);
        scanf("%lf", &results[i]);
        matrix[i][N] = results[i]; 
    }

    for (i = 0; i < N; i++) {
        double max = fabs(matrix[i][i]);
        int maxRow = i;
        for (k = i + 1; k < N; k++) {
            if (fabs(matrix[k][i]) > max) {
                max = fabs(matrix[k][i]);
                maxRow = k;
            }
        }

        for (k = i; k < N + 1; k++) {
            double temp = matrix[maxRow][k];
            matrix[maxRow][k] = matrix[i][k];
            matrix[i][k] = temp;
        }

        for (k = i + 1; k < N; k++) {
            double c = -matrix[k][i] / matrix[i][i];
            for (j = i; j < N + 1; j++) {
                if (i == j) {
                    matrix[k][j] = 0;
                } else {
                    matrix[k][j] += c * matrix[i][j];
                }
            }
        }
    }

    for (i = N - 1; i >= 0; i--) {
        solutions[i] = matrix[i][N] / matrix[i][i];
        for (k = i - 1; k >= 0; k--) {
            matrix[k][N] -= matrix[k][i] * solutions[i];
        }
    }

    printf("Çözümler:\n");
    for (i = 0; i < N; i++) {
        printf("x%d = %lf\n", i + 1, solutions[i]);
    }

    for (i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(results);
    free(solutions);
}

void inverseMatrix()
{
    int N;
    int i, j, k;

    printf("Matris için N sayısını girin: ");
    scanf("%d", &N);

    double **matrix = (double **)malloc(N * sizeof(double *));
    double **inverse = (double **)malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
    {
        matrix[i] = (double *)malloc(N * sizeof(double));
        inverse[i] = (double *)malloc(N * sizeof(double));
    }

    printf("Matris elemanlarını girin:\n");
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            printf("Matris[%d][%d]: ", i, j);
            scanf("%lf", &matrix[i][j]);
            if (i == j)
                inverse[i][j] = 1.0;
            else
                inverse[i][j] = 0.0;
        }
    }

    // Gauss-Jordan eliminasyonu
    for (i = 0; i < N; i++)
    {
        double temp = matrix[i][i];
        for (j = 0; j < N; j++)
        {
            matrix[i][j] /= temp;
            inverse[i][j] /= temp;
        }
        for (j = 0; j < N; j++)
        {
            if (i != j)
            {
                temp = matrix[j][i];
                for (k = 0; k < N; k++)
                {
                    matrix[j][k] -= matrix[i][k] * temp;
                    inverse[j][k] -= inverse[i][k] * temp;
                }
            }
        }
    }

    printf("Orijinal matris:\n");
    printMatrix(matrix, N, N);
    printf("Ters matris:\n");
    printMatrix(inverse, N, N);

    for (i = 0; i < N; i++)
    {
        free(matrix[i]);
        free(inverse[i]);
    }
    free(matrix);
    free(inverse);
}

void printMatrix(double **matrix, int rows, int cols) {
    int i, j;
    printf("matris:\n");
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

void trapezoidalMethod()
{
    double a, b, h, sum, counter, lastRoot, n, currentX;
    int i;

    char *expression = getExpression();
    printExpression(expression);
    printf("Starting point of the interval (a):");
    scanf("%lf", &a);
    printf("Ending point of the interval (b):");
    scanf("%lf", &b);
    printf("Enter the amount of subdivisions(n):");
    scanf("%lf", &n);

    sum = (calcExpression(expression, a) + calcExpression(expression, b)) / 2;
    currentX = a;
    h = (b - a) / n;
    for (i = 0; i < n - 1; i++)
    {
        currentX += h;
        sum = sum + (calcExpression(expression, currentX));
    }
    printf("Calculated integral for trapezoidal method:");
    printf("%lf", sum * h);
}

void simpsonMethod()
{
    double a, b, h, currentX, sum, counter, lastRoot, n;
    int i;

    char *expression = getExpression();
    printExpression(expression);
    printf("Starting point of the interval (a):");
    scanf("%lf", &a);
    printf("Ending point of the interval (b):");
    scanf("%lf", &b);
    printf("Enter the amount of subdivisions(n) (For both simpson 1/3 and 3/8 to converge, the amount of subdivisions has to be a factor of 6):");
    scanf("%lf", &n);

    // simpson 1/3
    sum = calcExpression(expression, a) + calcExpression(expression, b);
    lastRoot = a;
    h = (b - a) / n;
    for (i = 0; i < n - 1; i++)
    {
        currentX = lastRoot + h;
        if (i % 2 == 0)
        {
            sum = sum + 4 * (calcExpression(expression, currentX));
        }
        else
        {
            sum = sum + 2 * (calcExpression(expression, currentX));
        }
        lastRoot = currentX;
    }
    printf("Calculated integral for simpson 1/3:");
    printf("%lf", sum * h / 3);

    // simpson 3/8
    sum = calcExpression(expression, a) + calcExpression(expression, b);
    lastRoot = a;
    h = (b - a) / n;
    for (i = 0; i < n - 1; i++)
    {
        currentX = lastRoot + h;
        if (i % 3 == 2)
        {
            sum = sum + 2 * (calcExpression(expression, currentX));
        }
        else
        {
            sum = sum + 3 * (calcExpression(expression, currentX));
        }
        lastRoot = currentX;
    }
    printf("\nCalculated integral for simpson 3/8:");
    printf("%lf", sum * 3 * h / 8);
}

double differentiation(char *expression, int differentiationType, double h, double x)
{
    double temp1, temp2, temp3, result = 0.0;
    if (differentiationType == 1)
    {
        temp1 = calcExpression(expression, x);
        temp2 = calcExpression(expression, x - h);
        temp3 = calcExpression(expression, x - (2 * h));
        result = ((-3 * temp1) + (4 * temp2) - temp3) / (2 * h);
    }
    else if (differentiationType == 2)
    {
        temp1 = calcExpression(expression, x + h);
        temp2 = calcExpression(expression, x - h);
        result = (temp1 - temp2) / (2 * h);
    }
    else if (differentiationType == 3)
    {
        temp1 = calcExpression(expression, x);
        temp2 = calcExpression(expression, x + h);
        temp3 = calcExpression(expression, x + (2 * h));
        result = ((-3 * temp1) + (4 * temp2) - temp3) / (2 * h);
    }

    return result;
}

void numericDerivative()
{
    double h, x;
    int i, j, derivativeType;

    char *expression = getExpression();
    printExpression(expression);
    printf("1-Backward differentiation\n2-nCentral differentiation\n3-Forwward differentiation\n");
    scanf("%d", &derivativeType);
    printf("h: ");
    scanf("%lf", &h);
    printf("x: ");
    scanf("%lf", &x);

    double result = differentiation(expression, derivativeType, h, x);
    printf("result: %lf\n", result);
}

void newtonRaphsonMethod()
{
    double start, x, error, expectedError, realValue, result, tmp, tmp2;
    int i, equationCount, iterMax, iter = 0, errorType;

    // printf("\nenter element count: ");
    // scanf("%d", &equationCount);

    // POLY *equation = getEquation(equationCount);
    // printEquation(equation, equationCount);
    // POLY *derivative = getDerivative(equation, equationCount);
    // printEquation(derivative, equationCount);

    char *expression = getExpression();
    printExpression(expression);

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
        tmp = calcExpression(expression, x) / differentiation(expression, 2, 0.01, x);
        tmp2 = x;
        x = x - tmp;
        result = calcExpression(expression, x);
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

    free(expression);
}

void regulaFalsiMethod()
{
    double start, left, right, mid, error, expectedError, realValue, result, tmp;
    int i, equationCount, iterMax, iter = 0, errorType;

    // printf("\nenter element count: ");
    // scanf("%d", &equationCount);

    // POLY *equation = getEquation(equationCount);
    // printEquation(equation, equationCount);

    char *expression = getExpression();
    printExpression(expression);

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

    if (calcExpression(expression, left) * calcExpression(expression, right) > 0)
    {
        printf("invalid range ");
    }
    else
    {
        do
        {
            iter++;
            // PROBLEM
            tmp = calcExpression(expression, right) - calcExpression(expression, left);
            mid = left - (calcExpression(expression, left) * ((right - left) / tmp));
            result = calcExpression(expression, mid);
            // printf("tmp: %lf  mid: %lf  result: %lf\n", tmp, mid, result);
            if (calcExpression(expression, mid) * calcExpression(expression, left) < 0)
            {
                right = mid;
            }
            else if (calcExpression(expression, mid) * calcExpression(expression, right) < 0)
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

    free(expression);
}

void bisectionMethod()
{
    double start, left, right, mid, error, expectedError, realValue, result;
    int i, equationCount, iterMax, iter = 0, errorType;

    // printf("\nenter element count: ");
    // scanf("%d", &equationCount);

    // POLY *equation = getEquation(equationCount);
    // printEquation(equation, equationCount);

    char *expression = getExpression();
    printExpression(expression);

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

    if (calcExpression(expression, left) * calcExpression(expression, right) > 0)
    {
        printf("invalid range ");
    }
    else
    {
        do
        {
            iter++;
            mid = (left + right) / 2;
            result = calcExpression(expression, mid);
            if (calcExpression(expression, mid) * calcExpression(expression, left) < 0)
            {
                right = mid;
            }
            else if (calcExpression(expression, mid) * calcExpression(expression, right) < 0)
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

    free(expression);
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

void printExpression(char *expression)
{
    printf("'%s'\n", expression);
}

double calcExpression(char *expression, double x)
{
    int index = 0;
    double result = parse_expression(expression, &index, x);

    return result;
}

char *getExpression()
{
    char *expression = (char *)malloc(256 * sizeof(char));
    printf("Eşitsizlik gir: ");
    scanf("%s", expression);

    return expression;
}

double parse_exponent(char *expr, int *index, double x)
{
    double base = parse_factor(expr, index, x);

    while (expr[*index] == '^')
    {
        (*index)++; // Skip '^'
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

double parse_factor(char *expr, int *index, double x)
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

double parse_term(char *expr, int *index, double x)
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

double parse_expression(char *expr, int *index, double x)
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
    printf("Gauss elimination:5\n");
    printf("Gauss seidel Method:6\n");
    printf("Numerical Differentiation:7\n");
    printf("Simpson Method:8\n");
    printf("Trapezoidal Method:9\n");
    printf("Gregory Newton Interpolation:10\n");
    printf("Method of choice: ");
    scanf("%d", &method);
    return method;
}