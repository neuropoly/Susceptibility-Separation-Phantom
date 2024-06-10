
var clone =  function(x) {return JSON.parse(JSON.stringify(x)) }

var kmath = {}


module.exports = kmath;


kmath.matrix = function(v)
{
	var data;
	var sizes;
	if (typeof v[0] == 'number')
	{
		data =  v.slice(0);
		sizes = [v.length];
	}
	else
	{
		if (typeof v == 'object' && v._data != undefined)
		{
	 	    v = clone(v);
			return v;
		}

		data = [];
		sizes = [v.length, v[0].length];
		for (var k = 0 ; k < v.length ;k++)
			data[k] = v[k].slice(0);
	}
	return {_data:data,_size:sizes};
}

kmath.Transpose = function(v)
{
	if (v._data == undefined)
		v = kmath.matrix(v);
	else
		v = clone(v);

	if (v._size.length == 1)
	{
		v._size = [1,v._size[0]];
		v._data = [v._data];
	}

	var n = {_data:[],_size:[v._size[1],v._size[0]]}
	for (var k = 0; k< n._size[0];k++)
	{
		n._data[k] = [];
		for (var j = 0; j< n._size[1];j++)
		{
			n._data[k][j] = v._data[j][k];
		}
	}
    return n;

}
kmath.transpose = kmath.Transpose;


kmath.diag = function(v)
{

	var n = v.length;
	var data = [];
	for (var i = 0; i < n; i++)
	{
		data[i] = [];
		for (var j = 0; j < n; j++)
			data[i][j] = 0;
	}
	for (var i = 0; i < n; i++)
		data[i][i] = v[i];

	return kmath.matrix(data);
}

kmath.inv = function(M)
{
	M = kmath.matrix(M);
	var E = kmath.eye(M._size[0]);
	for (var k = 0; k < E._size[0];k++)
	{
		E._data[k] = kmath.mdiv(M,kmath.matrix(E._data[k]))._data;
	}

    return kmath.Transpose(E);
}


kmath.eye = function(n)
{
	var data = [];
	for (var i = 0; i < n; i++)
		data[i] = 1;
	return kmath.diag(data)
}

kmath.pow = Math.pow;
kmath.abs = Math.abs;
kmath.sqrt = Math.sqrt;
kmath.round = function(x)
{
	if (Array.isArray(x))
	{
		x = x.slice(0);
		for (var k = 0; k < x.length;k++)
		{
			x[k] = Math.round(x[k]);
		}
		return x;
	}
	else
		return Math.round(x);

}
kmath.floor = Math.floor;
kmath.ceil = Math.ceil;
kmath.exp = Math.exp;


kmath.sign = function(x) {
	return (x>=0)?1:-1; }


kmath.sum = function(x)
{
	var s = 0;
	for (var k = 0;k < x.length;k++)
	{
		s += x[k];
	}
	return s;

}


kmath.maxEV = function(M)
{
	var n = M._data.length;
	var v = new Array(n).fill(0); v[0] = 1;
	v = kmath.matrix(v);
	var old = v;
	var ev;
	for (var k = 0;k < 20; k++)
	{
		v = kmath.multiply(M,v);
		ev = kmath.normalize(v);
		var dif = 0;
		for (var j = 0; j < n;j++)
			dif += (v._data[j]-old._data[j])*(v._data[j]-old._data[j]);
		old = v;
		if (dif < 0.00000001)
			break;
	}

	return {v:v,ev:ev}	;

}

kmath.normalize = function (v)
{
	var n = 0;
	for (var k = 0; k < v._data.length;k++)
		n += v._data[k]*v._data[k];
	n = kmath.sqrt(n);
	for (var k = 0; k < v._data.length;k++)
		v._data[k] /= n;
	return n;
}


kmath.norm = function (v)
{
	var n = 0;
	for (var k = 0; k < v.length;k++)
	{
		n += v[k]*v[k];
	}
	n = kmath.sqrt(n);
	return n;
}



kmath.mdiv = function(M,b)
{
	var abs = Math.abs;
    M = clone(M);

    var A = M._data;
    var y = b._data.slice(0);

    var r = clone(b);
    r._data = gauss(A,y);

    return r;


	function array_fill(i, n, v) {
		var a = [];
		for (; i < n; i++) {
			a.push(v);
		}
		return a;
	}

	function gauss(A, x) {

		var i, k, j;

		// Just make a single matrix
		for (i=0; i < A.length; i++) {
			A[i].push(x[i]);
		}
		var n = A.length;

		for (i=0; i < n; i++) {
			// Search for maximum in this column
			var maxEl = abs(A[i][i]),
				maxRow = i;
			for (k=i+1; k < n; k++) {
				if (abs(A[k][i]) > maxEl) {
					maxEl = abs(A[k][i]);
					maxRow = k;
				}
			}


			// Swap maximum row with current row (column by column)
			for (k=i; k < n+1; k++) {
				var tmp = A[maxRow][k];
				A[maxRow][k] = A[i][k];
				A[i][k] = tmp;
			}

			// Make all rows below this one 0 in current column
			for (k=i+1; k < n; k++) {
				var c = -A[k][i]/A[i][i];
				for (j=i; j < n+1; j++) {
					if (i===j) {
						A[k][j] = 0;
					} else {
						A[k][j] += c * A[i][j];
					}
				}
			}
		}

		// Solve equation Ax=b for an upper triangular matrix A
		var z = array_fill(0, n, 0);
		for (i=n-1; i > -1; i--) {
			z[i] = A[i][n]/A[i][i];
			for (k=i-1; k > -1; k--) {
				A[k][n] -= A[k][i] * z[i];
			}
		}

		return z;
	}
}

kmath.multiply = function(a,b)
{
	if (Array.isArray(a))
		a = kmath.matrix(a);
	if (Array.isArray(b))
		b = kmath.matrix(b);

	if (typeof a == 'number')
	{
		var r = kmath.matrix(b._data);
		if (r._size.length == 1)
			for (var k = 0; k < r._size[0];k++)
				r._data[k] = r._data[k] * a;
		else
		{
			for (var k = 0; k < r._size[0];k++)
				for (var j = 0; j < r._size[1];j++)
					r._data[k][j] = r._data[k][j] * a;
		}
		return r;
	}
	else if (typeof b == 'number')
		return kmath.multiply(b,a);
	else
	{
		var data =  [];
		var trans = false;
		if (b._size.length == 1)
		{
			b = kmath.Transpose(b);
			trans = true;
		}

		for (var k=0;k < a._size[0];k++)
		{
			data[k] = [];
			for (var j=0;j < b._size[1];j++)
			{
				var v = 0;
				for (var i = 0; i < a._size[1];i++)
					v += a._data[k][i]*b._data[i][j];
				data[k][j] = v;
			}
		}

		var r =  {_data:data,_size:[a._size[0],b._size[1]]};
		if (trans)
		{
			r = kmath.Transpose(r);
			r._data = r._data[0];
			r._size = [r._size[1]];
			return r;
		}
	    else
	    	return r;


	}


}

kmath.max = function(a,b)
{
	if (b==undefined)
	{
		if (a._data != undefined)
			a = a._data

		var m = a[0]
		for (var k = 1; k < a.length; k++)
		{
			if (m<a[k])
				m = a[k];
		}
		return m;
	}
	else
	{
		if (a>b)
			return a;
		else
			return b;
	}
}


kmath.min = function(a,b)
{
	if (b==undefined)
	{
		if (a._data != undefined)
			a = a._data

		var m = a[0]
		for (var k = 1; k < a.length; k++)
		{
			if (m>a[k])
				m = a[k];
		}
		return m;
	}
	else
	{
		if (a<b)
			return a;
		else
			return b;
	}
}


kmath.add = function(a,b)
{
	if (typeof a == 'number')
		return kmath.add(b,a);

	a = kmath.matrix(a);

	if (typeof b == 'number')
	{
		if (a._size.length == 1)
		{
			for (var k = 0; k < a._size[0];k++)
				a._data[k] += b;
		}
		else
		{
			for (var k = 0; k < a._size[0];k++)
				for (var j = 0; j < a._size[1]; j++)
					a._data[k][j] += b;
		}
		return a;
	}
	else
	{
		if (a._size.length == 1)
		{
			for (var k = 0; k < a._size[0];k++)
				a._data[k] += b._data[k];
		}
		else
		{
			for (var k = 0; k < a._size[0];k++)
				for (var j = 0; j < a._size[1]; j++)
					a._data[k][j] += b._data[k][j];
		}

		return a;
	}



}



kmath.det = function(M)
{
	M = kmath.matrix(M);
	M = kmath.Transpose(M);
    var n = M._size[0];
	var a= M._data.reduce(function(a, b) {
	  return a.concat(b);
	});

	return Det();


	function Det()
	{
		if (n==1) return a[0];
		if (n==2) return a[0]*a[3]-a[1]*a[2];
		if (n==3) return a[0]*(a[4]*a[8] - a[5]*a[7])+a[1]*(a[5]*a[6]-a[3]*a[8])+a[2]*(a[3]*a[7]-a[4]*a[6]);
		if (n==4)
		{
			d= a[0]*(a[5]*(a[10]*a[15]-a[11]*a[14])+a[6]*(a[11]*a[13]-a[9]*a[15])+a[7]*(a[9]*a[14]-a[10]*a[13]));
			d-=a[1]*(a[4]*(a[10]*a[15]-a[11]*a[14])+a[6]*(a[11]*a[12]-a[8]*a[15])+a[7]*(a[8]*a[14]-a[10]*a[12]));
			d+=a[2]*(a[4]*(a[9]*a[15]-a[11]*a[13])+a[5]*(a[11]*a[12]-a[8]*a[15])+a[7]*(a[8]*a[13]-a[9]*a[12]))-a[3]*(a[4]*(a[9]*a[14]-a[10]*a[13])+a[5]*(a[10]*a[12]-a[8]*a[14])+a[6]*(a[8]*a[13]-a[9]*a[12]));
			return d;
		}
		detGLSL(n);
		var p=1.0;
		for (i=0; i<n; i++)
		{
			p*=a[i*(n+1)];
		}
		return p;
	}

}

math = kmath;
