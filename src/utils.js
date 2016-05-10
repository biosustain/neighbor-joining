import * as TimSort from 'timsort';

export function allocateSquareMatrix(n, value=null) {
    let a = new Array(n);
    for (let i = 0; i < n; i++) {
        a[i] = new Array(n);
        if(value !== null) a[i].fill(value);
    }
    return a;
}

export function arrayCopy(a) {
    let b = new Array(a.length),
        i = a.length;
    while(i--) { b[i] = a[i]; }
    return b;
}

export function sumRows(a) {
    let sum,
        n = a.length,
        sums = new Array(n);

    for (let i = 0; i < n; i++) {
        sum = 0;
        for (let j = 0; j < n; j++) {
            if (a[i][j] === undefined) continue;
            sum += a[i][j];
        }
        sums[i] = sum;
    }

    return sums;
}

export function sortWithIndices(toSort, skip=-1, timsort=false) {
    var n = toSort.length;
    var indexCopy = new Array(n);
    var valueCopy = new Array(n);
    var i2 = 0;

    for (var i = 0; i < n; i++) {
        if (toSort[i] === -1 || i === skip) continue;
        indexCopy[i2] = i;
        valueCopy[i2++] = toSort[i];
    }
    indexCopy.length = i2;
    valueCopy.length = i2;

    if (timsort) {
        TimSort.sort(indexCopy, (a, b) => toSort[a] - toSort[b]);
    }
    else {
        indexCopy.sort((a, b) => toSort[a] - toSort[b]);
    }

    TimSort.sort(indexCopy,function(left, right) {
        return toSort[left] - toSort[right];
    });

    valueCopy.sortIndices = indexCopy;
    for (var j = 0; j < i2; j++) {
        valueCopy[j] = toSort[indexCopy[j]];
    }
    return valueCopy;
}
