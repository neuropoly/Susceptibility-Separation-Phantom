//     unring.js
//
//     Elias Kellner, 
//     Marco Reisert,
//     Paul Reggentin
//     Medical Physics, University Medical Center Freiburg

/*  
    LICENCE
        
    The Software remains the property of the University Medical Center 
    Freiburg ("the University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. */





if (typeof module != "undefined")
{

    var math    = require('./kmath.js');
    var fs      = require("fs");
    var nifti   = require('./nifti.js');
    //var FFT = require("fft");
    var KissFFT = require("./kissFFT_main.js");

//     module.exports =
//     {
//         unringer:undefined
//     }
}



/////////////////////////////////////////////////////
//  Command line call
/////////////////////////////////////////////////////
if(typeof process != "undefined")
{

    //  use --key values pairs
    //  for params see  the "unring_nifti" function

    var params = {};

    var args = process.argv.slice(2);
    if(args.length > 0)
    {
        for(var k=0; k < args.length; k+=2)
        {
            if(args[k].substr(0,2) !="--")
            {
               console.log("Error parsing input `"+ args[k] +"`. Must be given as [--key value] pair")
               
                
			   //process.exit(1);
			   //return;
            }
            else
            {
               if(args[k] =="--help")
               {
                   printhelp();
                   //return;
               } 
               params[args[k].substr(2)] = parseInt(args[k+1]); 
            }

        }
    }
    else
    {
        printhelp();
    }

    function printhelp()
    {
        var nl = "\n\t\t\t ";
        var msg = "";
        msg += "Usage: nodejs unring.js";
        msg += "  --in <file_in>";
        msg += (nl + "--out <file_in> (default: <in>_unring)");
        msg += (nl + "--nshifts x (default: 20)");
        msg += (nl + "--minW    x (default: 1)");
        msg += (nl + "--maxW    x (default: 3)");

        msg += (nl + "  ===== test / debug  params =====");
        msg += (nl + "--maxnumslices     x (default: -1)");
        msg += (nl + "--maxnumtimepoints x (default: -1)");
        msg += (nl + "--interleaved      x (default:  0)");
        console.log(msg)
        return;
    }
    
    if(params.in == undefined)
        params.in  = 't2zoom_test.nii';
        
    if(params.out == undefined)
        params.out = params.in.substr(0, params.in.length-3) + "_unring.nii";

    console.log('Unringing file ... `' + params.in +"`");

    unring_file(params, function(){} );
    
    console.log('Finished. Results written to `' + params.out + '`');

}




//////////////////////////////////////////////////
// command line wrapper (nodejs)
//////////////////////////////////////////////////
function unring_file(params)
{
    var raw = fs.readFileSync(params.in);
    raw.buffer = raw.buffer || raw;
    var nii = parse(raw.buffer);
    
    unring_nifti(nii, params)

    fs.writeFileSync(params.out, raw);

}



//////////////////////////////////////////////////
// main function
//////////////////////////////////////////////////
function unring_nifti(nii, params_in)
{
    var params = {};
    
   
    // ********** main params
    params.nshifts = 20;
    params.minW = 1;
    params.maxW = 2;

    // ********** debug / test params
    params.maxnumslices = -1;
    params.maxnumtimepoints = -1;
    params.interleaved  = false;

    params.whendone = function(){};
    params.onprogress = console.log;
   


    
    // *********** extend with the input params
    if(params_in !=undefined)
    {
        for(var p in params)
        {
            if(params_in[p] != undefined)
                if( typeof(params_in[p]) != "function" && typeof params_in[p] != "boolean" ) 
                    params[p] = parseInt( params_in[p] ); 
                else
                     params[p] = ( params_in[p] ); 
        }
    }


    // nifti buffer will be overwritten, so give this function a copy of a nifti, if desired

    var data = nii.data;
    var sz = nii.sizes;

    var raw_slice    = new Float64Array(sz[0]*sz[1]);
    var unrung_slice = new Float64Array(sz[0]*sz[1]*2);
    var unrung_int16 = new Int16Array(sz[0]*sz[1]);
    
    if(sz[3] == undefined)
        sz[3] = 1;
    
    var maxtp = params.maxnumtimepoints==-1?sz[3]:params.maxnumtimepoints;
    if(maxtp > sz[3])
        maxtp = sz[3];
    
    var maxslices  = params.maxnumslices==-1?sz[2]:params.maxnumslices;
    if(maxslices > sz[2])
        maxslices = sz[2];
    
    // time measure nodejs /javascript
    if(typeof performance != "undefined")
    {
        var gettime = function(){ return performance.now()};
    }
    else
    {
         var gettime = function(){ return new Date().getTime()};
    }
    
    console.log(" : Parameters: nshifts " +  params.nshifts + ", maxW " + params.minW + ", maxW " + params.maxW + ", maxnumslices " + maxslices +  ", maxnumtimepoints " + maxtp );

    //test / debug: unring only each second slice and interleave with original
    if(params.interleaved)
    {
        //**** outermost loop ist timepoint (b-vals for diffusion data)
        for( var timepoint=0; timepoint < maxtp; timepoint++)
        {
            for(var slice = 1; slice < maxslices; slice +=2 )
            {
                
                var st = gettime();
                var offsetbefore = (slice -1)*sz[0]*sz[1] + timepoint*sz[0]*sz[1]*sz[2];
                var offset = (slice)*sz[0]*sz[1] + timepoint*sz[0]*sz[1]*sz[2];

                raw_slice = nii.data.slice(offsetbefore, offset + sz[0]*sz[1]);
                unring2d(raw_slice,false,unrung_slice,sz, params.nshifts, params.minW, params.maxW);

                //convert unrung data to int16 for output
                for(var i = 0; i < unrung_int16.length; i++)
                {
                    nii.data[i+offset] = unrung_slice[2*i] + .5;
                }

                var en = gettime();
                
                var msg ="";
                if(sz[3] > 1)
                    msg += "Unringing volume " +timepoint  + "/" +sz[3] + " | ";
                msg += "slice " +slice  + "/" +sz[2] + ' in ' +  Math.round(en-st) + 'ms.' ;
                params.onprogress( msg );

            }//for slice
        }// for time

    }
    else // normal mode
    {    
        //**** outermost loop ist timepoint (b-vals for diffusion data)
        for( var timepoint=0; timepoint < maxtp; timepoint++)
        {
            for(var slice = 0; slice < maxslices; slice++)
            {
                var st = gettime();
                var offset = slice*sz[0]*sz[1] + timepoint*sz[0]*sz[1]*sz[2];

                //selslice = 4* sz[0]*sz[1];

                raw_slice = nii.data.slice(offset, offset + sz[0]*sz[1]);

                unring2d(raw_slice,false,unrung_slice,sz, params.nshifts, params.minW, params.maxW);

                //convert unrung data to int16 for output
                for(var i = 0; i < unrung_int16.length; i++)
                {
                    nii.data[i+offset] = unrung_slice[2*i] + .5;
                }

                 var en = gettime();

                var msg ="";
                if(sz[3] > 1)
                    msg += "Unringing volume " +timepoint  + "/" +sz[3] + " | ";
                msg += "slice " +slice  + "/" +sz[2] + ' in ' +  Math.round(en-st) + 'ms.' ;
                params.onprogress( msg );

            }//for slice
        }// for time
    }

        
    // final callback
    params.whendone();
    console.log(" Unringing done." );

    return nii;
}



/*
    img_in: 2D image, with possible ringing, stored as 1D Float64Array 
    img_out: 2D image stored as 1D Float64Array which will store the unrung image
    dim: 2-elt array with the dimensions of the input image
    nshifts: number of shifts to try
    minW: minimum of the shift window
    maxW: maximum of the shift window
*/
function unring2d(img_in, is_cplx, img_out, dim, nshifts, minW,maxW)
{

    //create the fourier transforms of the image and its transpose
    var fac = 2;
    if(is_cplx)
    {
        fac = 1;
    }
    var img_in_ft = new Float64Array(img_in.length*fac);
    var img_in_t_ft = new Float64Array(img_in.length*fac);

    fft_2d(img_in,dim,img_in_ft,is_cplx);
    transpose_cplx(img_in_ft,dim,img_in_t_ft);
    //(the 2d fft of the transpose of X i the hermitian transpose of the 2d fft of x)

    //create and apply saddle filter
    for(var j = 0; j < dim[1]; j++)
    {
        var cj = (1+Math.cos(2*Math.PI*j/dim[1]))*0.5;
        for(var i = 0; i < dim[0]; i++){
            var ci = (1+Math.cos(2*Math.PI*i/dim[0]))*0.5;           
            var eps = .000000000000001
            var scale = 1/(ci+cj+eps);//the denominator of the filter value, eps to avoid /0 error
            img_in_ft[2*(i+j*dim[0])] = img_in_ft[2*(i+j*dim[0])] * cj * scale;
            img_in_ft[2*(i+j*dim[0]) + 1] = img_in_ft[2*(i+j*dim[0]) + 1] * cj * scale;
            img_in_t_ft[2*(j + i*dim[1])] = img_in_t_ft[2*(j + i*dim[1])] * ci * scale;  
            img_in_t_ft[2*(j + i*dim[1]) + 1] = img_in_t_ft[2*(j + i*dim[1]) + 1] * ci * scale;  
        }//for i
    }//for j

    //convert the filtered images back to spatial domain
    var filt1 = new Float64Array(img_in.length*2);
    var filt2 = new Float64Array(img_in.length*2);
    ifft_2d(img_in_ft,dim,filt1);
    ifft_2d(img_in_t_ft,dim_t,filt2);

    //unring the filtered images
    unring1d(filt1,dim[0],dim[1],nshifts,minW,maxW,function(){});
    unring1d(filt2,dim[1],dim[0],nshifts,minW,maxW,function(){});

    //add the two filtered images together and store them in img_out
    for(var i = 0; i < dim[0]; i++)
    {
        for(var j = 0; j < dim[1]; j++){
            img_out[2*(i+j*dim[0])] = filt1[2*(i+j*dim[0])] + filt2[2*(j + i*dim[1])];
            img_out[2*(i+j*dim[0]) + 1] = filt1[2*(i+j*dim[0]) + 1] + filt2[2*(j + i*dim[1]) + 1];
        }
    }
}//unring2d

/*
    transpose a real-valued 2-d image 
    img: 1-D vector storing the original image
    dim: the dimensions of the 2-d image stored in img vector
    img_t: 1-D vector that will hold the transposed image 
*/
function transpose_re(img,dim,img_t)
{
    for(var i = 0; i < dim[0]; i++){
        for(var j = 0; j < dim[1]; j++){
            img_t[j+i*dim[1]] = img[i+j*dim[0]];
        }
    }
}

/*
    transpose a real-valued 2-d image 
    img: 1-D vector storing the original image 
            in the form [real0 imag0 real1 imag1 ... ]
    dim: the dimensions of the 2-d image stored in img vector
    img_t: 1-D vector that will hold the transposed image 
*/
function transpose_cplx(img,dim,img_t)
{
    for(var i = 0; i < dim[0]; i++)
    {
        for(var j = 0; j < dim[1]; j++)
        {
            img_t[2*(j+i*dim[1])] = img[2*(i+j*dim[0])];
            img_t[2*(j+i*dim[1]) + 1] = img[2*(i+j*dim[0]) + 1];
        }
    }
}

/*
    insert a source vector into a destination vector
    src: the vector to be inserted
    dest: the vector that src is inserted into
    n: number of elements of src to insert
    offset: the position in dest that src is inserted into
*/
function copy_into(src,dest,n,offset)
{
    for(var i = 0; i < n; i++)
    {
        dest[offset+i] = src[i];
    }
}

/*
    perform a 2-dimensional fourier transform 
    img: the image to be transformed, stored as a 1D vector. 
        can be real or complex, stored as [real0 imag0 real1 imag1 ...]
    dim: the dimensions of the 2D image that 'img' represents
    img_ft: the vector that will hold the fourier transform of img
    is_cplx: boolean flag for type of img
*/
function fft_2d(img, dim, img_ft, is_cplx)
{
    //create fft objects
    //fft0 = new FFT.complex(dim[0],false);
    //fft1 = new FFT.complex(dim[1],false);

    var kissFFT0 = new KissFFT.FFT(dim[0]);
    var kissFFT1 = new KissFFT.FFT(dim[1]);


    //create vectors that will store intermediate steps
    var frow = new Float64Array(2*dim[0]);
    var intermed = new Float64Array(2*dim[0]*dim[1]);
    var intermed_t = new Float64Array(2*dim[0]*dim[1]);
    var fft_trans = new Float64Array(2*dim[0]*dim[1]);


    // make image complex if not already, makes it easier, only complex fft
    if(!is_cplx)
    {
        var img2 = new Float64Array(dim[0]*dim[1]*2);
        
        // not working in all nodejs versions ... so zero fill manually
        //img2.fill(0);
        for(var k=0; k< img2.length; k++) 
            img2[k] = 0;

        for(var k=0; k< img.length; k++)
            img2[2*k] = img[k];
        img = img2;
        var fft_type = 'complex';
    }

    //perform 1D fft on each of the rows of img
    var row_len = 2*dim[0];
    for(var i = 0; i < dim[1]; i++)
    {
        //fft0.simple(frow,img.slice(i*row_len,(i+1)*row_len),fft_type);
        frow = kissFFT0.forward(img.slice(i*row_len,(i+1)*row_len));
        copy_into(frow,intermed,2*dim[0],i*2*dim[0]);
    }

    //transpose the row transformed array
    transpose_cplx(intermed,dim,intermed_t);

    //fft each of the rows (which are columns from the original image)
    var col_len = 2*dim[1];
    var fcol = new Float64Array(col_len);
    for(var i = 0; i < dim[0]; i++)
    {
        //fft1.simple(fcol,intermed_t.slice(i*col_len,(i+1)*col_len),'complex');
        fcol = kissFFT1.forward(intermed_t.slice(i*col_len,(i+1)*col_len));
        copy_into(fcol,fft_trans,col_len,i*col_len);
    }

    //transpose back to original dimentions
    dim_t = new Float64Array([dim[1], dim[0]]);
    transpose_cplx(fft_trans,dim_t,img_ft);
}

//inverse fourier transform
//same as fft_2d but all of the 1D ffts are inverse DFTs
function ifft_2d(kspace, dim, img_space)
{
//      fft0 = new FFT.complex(dim[0],true);
//      fft1 = new FFT.complex(dim[1],true);

    var ikissFFT0 = new KissFFT.FFT(dim[0]);
    var ikissFFT1 = new KissFFT.FFT(dim[1]);

    //fft each of the rows
    var frow = new Float64Array(2*dim[0]);
    var intermed = new Float64Array(2*dim[0]*dim[1]);
    var intermed_t = new Float64Array(2*dim[0]*dim[1]);
    var fft_trans = new Float64Array(2*dim[0]*dim[1]);

    var fft_type = 'complex';
    var row_len = 2*dim[0];

    for(var i = 0; i < dim[1]; i++)
    {
        //fft0.simple(frow,kspace.slice(i*row_len,(i+1)*row_len),fft_type);
        frow = ikissFFT0.inverse( kspace.slice(i*row_len,(i+1)*row_len)  );
        copy_into(frow,intermed,2*dim[0],i*2*dim[0]);
        
    }

    //transpose the row transformed
    transpose_cplx(intermed,dim,intermed_t);

    //fft each of the rows (which are columns from the original image)
    var col_len = 2*dim[1];
    var fcol = new Float64Array(col_len);
    for(var i = 0; i < dim[0]; i++)
    {
        //fft1.simple(fcol,intermed_t.slice(i*col_len,(i+1)*col_len),'complex');
        fcol = ikissFFT1.inverse(intermed_t.slice(i*col_len,(i+1)*col_len) )
        copy_into(fcol,fft_trans,col_len,i*col_len);
    }
    //transpose back to original dimentions
    dim_t = new Float64Array([dim[1], dim[0]]);
    transpose_cplx(fft_trans,dim_t,img_space);
    
    //perform scaling
    scale = 1/(dim[0]*dim[1]);
    for(var i = 0; i < img_space.length; i++)
    {
        img_space[i] = img_space[i]*scale; 
    }
}



/*
perform the unringing algorithm on each row of a 2-d array
ARGS: 
    data     - the numlines-by-n image array
    n        - number of pixels per row
    numlines - number of rows to analzye
    nshifts      - number of subvoxel shift values to try in one direction (total # shifts will be 2*nshifts+1)
    minW     - min of window
    maxW     - max of window
    callback - callback function
*/
function unring1d(data, n, numlines, nshifts, minW, maxW, callback)
{

   
    //create the length variables
    var totalShifts = 2*nshifts+1;
    var totalSamps = 2*n;
    
    //create the necessary variables
//     var fft  = new FFT.complex(n,false);
//     var ifft = new FFT.complex(n,true);

    var kissFFT = new KissFFT.FFT(totalSamps/2);
    var kissFFTi = new KissFFT.FFT(totalSamps/2);


    //this variable holds the frequencydomain version of a row, so that it can be shifted
    var rowFreq = new Float64Array(totalSamps);
    var rowFreqShift = new Float64Array(totalSamps);
    var rowTimeShift = new Float64Array(totalSamps);

    //these variables hold the rows shifted in time and freq domains
    var timeShifts = new Float64Array(totalSamps*totalShifts);
    //var freqShifts = new Float64Array(totalSamps*totalShifts);

    //these arrays hold the measure of ringing at successive shifts to the right and leftDiffs
    //ringing is measured as the sum of the pairwise differences of adjacent pixels
    var rightDiffs = new Float64Array(totalShifts);
    var leftDiffs = new Float64Array(totalShifts);


    //create the shift array
    var shifts = new Float64Array(totalShifts);
    shifts[0] = 0;
    for(var i = 0; i < nshifts; i++)
    {
        shifts[i+1] = (i+1);
        shifts[i+1+nshifts] = 0-(i+1);
    }
    
    //the angle corresponding to a certain linear shift
    var phi = new Float64Array(2);
    //the angle that increases for each successive frequencz
    var ang = new Float64Array(2);


    //var maxn;
    if(n%2 == 1) 
    {
        var maxn = (n-1)/2;
    } 
    else 
    {
        var maxn = n/2-1;
    }
    
    //generate the frequency ramp
    var freqRamps = new Float64Array(2*maxn*shifts.length);
    var ramplen = 2*maxn;
    for(s = 1; s < shifts.length; s++)
    {
        var p =  shifts[s] * Math.PI / (n*nshifts);  
        phi[0] = Math.cos(p);
        phi[1] = Math.sin(p);
        ang[0] = 1;
        ang[1] = 0;
        for(var i =0; i < maxn; i++)
        {
            var tmp = ang[0];
            ang[0] = phi[0]*ang[0] - phi[1]*ang[1];
            ang[1] = tmp*phi[1] + phi[0]*ang[1];
            freqRamps[s*ramplen + 2*i] = ang[0];
            freqRamps[s*ramplen + 2*i + 1] = ang[1];
        }
    }


    //apply shift in frequency domain, then in time domain
    for (var line = 0; line < numlines; line++)
    {
        var line_idx = line*totalSamps;

        var currentline = data.slice(line_idx,line_idx+totalSamps); 

         //DEBUG
        //start_time = performance.now();
        for(var i = 0; i < totalSamps; i++)
            timeShifts[i] = data[line_idx+i];
        
        

        //fft.simple(rowFreq,currentline,'complex');
        rowFreq =  kissFFT.forward( currentline );// / totalSamps;

        //shift the frequency data by each value in shifts
        //s is shift number
        //var angle_scale = Math.PI /(n*nshifts);
        // start from 1 -> leave first image untouched
        for(var s = 1; s < shifts.length; s++)
        {
            for(var j = 1; j < rowFreq.length; j++)
                rowFreqShift[j] = rowFreq[j];

            
            //set the dc term and nyquist frequency
            rowFreqShift[0] = rowFreq[0];
            rowFreqShift[1] = rowFreq[1];
            
            //if the Nquist frequency is included, set it to 0
            if(n%2==0) 
            {
                rowFreqShift[n] = 0;
                rowFreqShift[ n + 1] = 0;
            }
            
            var tmp;
            for(var i =0; i < maxn; i++)
            {
                var L = (i+1);
                var FR = 2*i;

                rowFreqShift[2*L] = freqRamps[s*ramplen + FR]*rowFreq[2*L] - freqRamps[s*ramplen + FR + 1]*rowFreq[2*L+1];
                rowFreqShift[2*L+1] = freqRamps[s*ramplen + FR]*rowFreq[2*L+1] + freqRamps[s*ramplen + FR + 1]*rowFreq[2*L];
                
                L = (n-1-i);
                rowFreqShift[2*L] = freqRamps[s*ramplen + FR]*rowFreq[2*L] + freqRamps[s*ramplen + FR + 1]*rowFreq[2*L+1];
                rowFreqShift[2*L+1] = freqRamps[s*ramplen + FR]*rowFreq[2*L+1] - freqRamps[s*ramplen + FR + 1]*rowFreq[2*L];
                
            }//for i (freq bin to shift)
            

            //perform the inverse fft to get signal shifted in time domain
            //ifft.simple(rowTimeShift,rowFreqShift,'complex');
            rowTimeShift = kissFFTi.inverse( rowFreqShift ) ;

            //copy the time shifted array into the dictionary of shifted values
            for(var i = 0; i < totalSamps; i++)
                timeShifts[s*totalSamps + i] = rowTimeShift[i] / n;
            
        }//for s (shift number)

        
        // calc the first pixel. For each following pixel, simple subtract tail and add new head. 
        // Only two operations instead of a loop
        for(var s = 0; s<totalShifts; ++s)
        {
            var offset = s*totalSamps;
            rightDiffs[s] = 0;
            leftDiffs[s] = 0;
            
            for(var d = minW; d <= maxW; d++)
            {
                var rightIdx1 = (2*d + totalSamps)  %totalSamps;
                var rightIdx2 = (rightIdx1+2)       %totalSamps;
                var leftIdx1  = (-2*d + totalSamps) %totalSamps;
                var leftIdx2  = (leftIdx1-2)        %totalSamps;
                
                var diff = timeShifts[offset + rightIdx1] - timeShifts[offset + rightIdx2];
                if (diff>0) 
                    rightDiffs[s] += diff;
                else 
                    rightDiffs[s] -=diff;
                

                diff = timeShifts[offset + leftIdx1] - timeShifts[offset + leftIdx2];
                if (diff>0) 
                    leftDiffs[s] += diff;
                else
                    leftDiffs[s] -= diff;

            }//for d

        }//for s

        
        var oldRight1, oldRight2, oldLeft1, oldLeft2, newRight1, newRight2, newLeft1, newLeft2, diff;
        //for each pixel find the minimal ringing measure & shift to that ringing measure
        for(var pix = 0; pix < n; pix++)
        {
            var minDiff = 99999999999;
            var minIndex = 0;

            for(var s = 0; s< totalShifts; s++)
            {

                var offset = s*totalSamps;        

                if (rightDiffs[s] < minDiff)
                {
                    minDiff = rightDiffs[s];
                    minIndex = s;
                } 
                if (leftDiffs[s] < minDiff)
                {
                    minDiff = leftDiffs[s];
                    minIndex = s;
                }

                //get the ringing measure for the successive pixel by removing the distance measure on one end and adding the distance measure on the other end.
                oldRight1 = (2*(pix+minW)+totalSamps)    % totalSamps;
                oldRight2 = (oldRight1 + 2)              % totalSamps;
                newRight1 = (2*(pix+maxW+1) + totalSamps)% totalSamps;
                newRight2 = (newRight1 + 2)              % totalSamps;

                oldLeft1 = (2*(pix - maxW) + totalSamps) %totalSamps;
                oldLeft2 = (oldLeft1 -2 + totalSamps)    %totalSamps;
                newLeft1 = (2*(pix-minW+1) + totalSamps) %totalSamps;
                newLeft2 = (newLeft1-2 + totalSamps )    %totalSamps;

                diff = timeShifts[offset + oldRight1] - timeShifts[offset + oldRight2];
                if(diff > 0)
                    rightDiffs[s] -= diff;
                else
                    rightDiffs[s] += diff;

                diff = timeShifts[offset + oldLeft1] - timeShifts[offset + oldLeft2];
                if(diff > 0)
                    leftDiffs[s] -= diff;
                else
                    leftDiffs[s] += diff;

                diff = timeShifts[offset + newRight1] - timeShifts[offset + newRight2];
                if(diff > 0)
                    rightDiffs[s] += diff;
                else
                    rightDiffs[s] -= diff;

                diff = timeShifts[offset + newLeft1]  - timeShifts[offset + newLeft2];
                if(diff > 0)
                    leftDiffs[s] += diff;
                else
                    leftDiffs[s] -= diff;

            }//for s
            
            var sh = shifts[minIndex]/(2*nshifts);
            var shift_offset = minIndex*totalSamps;
            if(sh>0)
            {
                var dif0re = timeShifts[shift_offset + (2*(pix - 1) + totalSamps)%totalSamps];
                var dif1re = timeShifts[shift_offset + (2*pix + totalSamps)%totalSamps];
                data[line_idx+2*pix] = (dif1re*(1-sh) + dif0re*sh) ;
            }
            else
            {
                var dif1re = timeShifts[shift_offset + (2*pix + totalSamps)%totalSamps];
                var dif2re = timeShifts[shift_offset + (2*(pix + 1) + totalSamps)%totalSamps];
                data[line_idx+2*pix] = (dif1re*(1+sh) - dif2re*sh) ;
            }
            // remove negative values, (we are at uint16 at the end)
            if(data[line_idx+2*pix] < 0)
            {
              data[line_idx+2*pix] = 0;
            }

         
        }//for pix
    }//for line

}

function exportnii(vec,filename)
{
    fs.writeFileSync(filename,vec.toString());
}
