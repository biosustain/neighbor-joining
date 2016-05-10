'use strict';

Object.defineProperty(exports, "__esModule", {
    value: true
});
exports.allocateSquareMatrix = allocateSquareMatrix;
exports.arrayCopy = arrayCopy;
exports.sumRows = sumRows;
exports.sortWithIndices = sortWithIndices;

var _timsort = require('timsort');

var TimSort = _interopRequireWildcard(_timsort);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function allocateSquareMatrix(n) {
    var value = arguments.length <= 1 || arguments[1] === undefined ? null : arguments[1];

    var a = new Array(n);
    for (var i = 0; i < n; i++) {
        a[i] = new Array(n);
        if (value !== null) a[i].fill(value);
    }
    return a;
}

function arrayCopy(a) {
    var b = new Array(a.length),
        i = a.length;
    while (i--) {
        b[i] = a[i];
    }
    return b;
}

function sumRows(a) {
    var sum = void 0,
        n = a.length,
        sums = new Array(n);

    for (var i = 0; i < n; i++) {
        sum = 0;
        for (var j = 0; j < n; j++) {
            if (a[i][j] === undefined) continue;
            sum += a[i][j];
        }
        sums[i] = sum;
    }

    return sums;
}

function sortWithIndices(toSort) {
    var skip = arguments.length <= 1 || arguments[1] === undefined ? -1 : arguments[1];
    var timsort = arguments.length <= 2 || arguments[2] === undefined ? false : arguments[2];

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
        TimSort.sort(indexCopy, function (a, b) {
            return toSort[a] - toSort[b];
        });
    } else {
        indexCopy.sort(function (a, b) {
            return toSort[a] - toSort[b];
        });
    }

    TimSort.sort(indexCopy, function (left, right) {
        return toSort[left] - toSort[right];
    });

    valueCopy.sortIndices = indexCopy;
    for (var j = 0; j < i2; j++) {
        valueCopy[j] = toSort[indexCopy[j]];
    }
    return valueCopy;
}