%Matrixs Operations
function result = matrixMultiplication(matrix1, matrix2)
    [m, n] = size(matrix1);
    [p, q] = size(matrix2);

    if n ~= p
        error("The number of columns in matrix1 should be equal to the number of rows in matrix2.");
    end

    result = zeros(m, q);

    for i = 1:m
        for j = 1:q
            for k = 1:n
                result(i, j) = result(i, j) + matrix1(i, k) * matrix2(k, j);
            end
        end
    end
end

% Example usage:
matrix1 = [1, 2, 3; 4, 5, 6];
matrix2 = [7, 8; 9, 10; 11, 12];

result = matrixMultiplication(matrix1, matrix2);
disp(result);
%Numerical Methods:
function root = bisectionMethod(f, a, b, epsilon, maxIterations)
    if sign(f(a)) == sign(f(b))
        error("The function values at a and b should have opposite signs.");
    end

    iteration = 0;

    while abs(b - a) > epsilon && iteration < maxIterations
        c = (a + b) / 2;

        if abs(f(c)) < epsilon
            root = c;
            return;
        end

        if sign(f(c)) == sign(f(a))
            a = c;
        else
            b = c;
        end

        iteration = iteration + 1;
    end

    root = (a + b) / 2;
end

% Example usage:
f = @(x) x^3 - x - 1;
a = 1;
b = 2;
epsilon = 1e-6;
maxIterations = 100;

root = bisectionMethod(f, a, b, epsilon, maxIterations);
disp(root);
%Sorting Algorithms:Bubble sort
function sortedArray = bubbleSort(array)
    n = length(array);

    for i = 1:n-1
        swapped = false;

        for j = 1:n-i
            if array(j) > array(j+1)
                temp = array(j);
                array(j) = array(j+1);
                array(j+1) = temp;
                swapped = true;
            end
        end

        if ~swapped
            break;
        end
    end

    sortedArray = array;
end

% Example usage:
arr = [7, 4, 9, 2, 1, 5];
sortedArr = bubbleSort(arr);
disp(sortedArr);
%Merge sort
function sortedArray = mergeSort(array)
    n = length(array);

    if n <= 1
        sortedArray = array;
        return;
    end

    mid = fix(n / 2);
    left = mergeSort(array(1:mid));
    right = mergeSort(array(mid+1:n));

    sortedArray = merge(left, right);
end

function mergedArray = merge(left, right)
    nLeft = length(left);
    nRight = length(right);
    mergedArray = zeros(1, nLeft + nRight);

    i = 1;
    j = 1;
    k = 1;

    while i <= nLeft && j <= nRight
        if left(i) <= right(j)
            mergedArray(k) = left(i);
            i = i + 1;
        else
            mergedArray(k) = right(j);
            j = j + 1;
        end
        k = k + 1;
    end

    while i <= nLeft
        mergedArray(k) = left(i);
        i = i + 1;
        k = k + 1;
    end

    while j <= nRight
        mergedArray(k) = right(j);
        j = j + 1;
        k = k + 1;
    end
end

% Example usage:
arr = [7, 4, 9, 2, 1, 5];
sortedArr = mergeSort(arr);
disp(sortedArr);
%function sortedArray = quickSort(array)
    n = length(array);

    if n <= 1
        sortedArray = array;
        return;
    end

    pivotIndex = fix(n / 2);
    pivot = array(pivotIndex);

    % Partition the array
    smaller = [];
    equal = [];
    larger = [];

    for i = 1:n
        if array(i) < pivot
            smaller = [smaller, array(i)];
        elseif array(i) == pivot
            equal = [equal, array(i)];
        else
            larger = [larger, array(i)];
        end
    end

    % Recursively sort the smaller and larger subarrays
    sortedSmaller = quickSort(smaller);
    sortedLarger = quickSort(larger);

    % Concatenate the sorted subarrays
    sortedArray = [sortedSmaller, equal, sortedLarger];
end

% Example usage:
arr = [7, 4, 9, 2, 1, 5];
sortedArr = quickSort(arr);
disp(sortedArr);
%Search
function index = binarySearch(array, target)
    n = length(array);
    left = 1;
    right = n;

    while left <= right
        mid = fix((left + right) / 2);

        if array(mid) == target
            index = mid;
            return;
        elseif array(mid) < target
            left = mid + 1;
        else
            right = mid - 1;
        end
    end

    % If the target is not found
    index = -1;
end

% Example usage:
arr = [1, 2, 4, 5, 7, 9];
target = 2;
index = binarySearch(arr, target);
disp(index);
% Recrussion
function result = factorial(n)
    if n == 0 || n == 1
        result = 1;
    else
        result = n * factorial(n - 1);
    end
end

% Example usage:
n = 5;
factorialResult = factorial(n);
disp(factorialResult);
%Palindrome
function isPalindrome = checkPalindrome(str)
    str = lower(str);
    n = length(str);

    % Remove non-alphanumeric characters
    str = regexprep(str, '[^a-z0-9]', '');

    % Check if the string is a palindrome
    isPalindrome = true;
    for i = 1:fix(n/2)
        if str(i) ~= str(n-i+1)
            isPalindrome = false;
            break;
        end
    end
end

% Example usage:
string1 = 'A man, a plan, a canal, Panama!';
isPalindrome1 = checkPalindrome(string1);
disp(isPalindrome1);

string2 = 'Hello, World!';
isPalindrome2 = checkPalindrome(string2);
disp(isPalindrome2);
% Mean, Median, Mode
function [meanValue, medianValue, modeValue] = calculateStatistics(data)
    % Calculate the mean
    meanValue = mean(data);

    % Calculate the median
    medianValue = median(data);

    % Calculate the mode
    modeValue = mode(data);
end

% Example usage:
dataset = [1, 2, 3, 4, 4, 5, 6, 6, 6, 7];
[meanResult, medianResult, modeResult] = calculateStatistics(dataset);

disp(meanResult);
disp(medianResult);
disp(modeResult);

