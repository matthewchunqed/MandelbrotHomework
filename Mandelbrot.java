/*
 * Note that this program uses the Java Vector API, which is only
 * available as an incubator feature in Java 16+. In order to compile
 * and run the program you need to use flag '--add-modules
 * jdk.incubator.vector'. For example:
 *
 * > javac --add-modules jdk.incubator.vector Mandelbrot.java
 *
 * In order to run any program that uses the Mandelbrot class, you
 * similarly have to add the same flag to the java command, e.g.:
 *
 * > java --add-modules jdk.incubator.vector MandelbrotTester
 */

import jdk.incubator.vector.*;
import java.util.Arrays;


// A class for computing the Mandelbrot set and escape times for
// points in the complex plane
public class Mandelbrot {

    // the maximum number of iterations of the system to consider 
    private float maxIter = 100.0F;

    // the squared distance to the origin for escape
    private float maxSquareModulus = 100.0F;

    // coordinates of the region to be computed
    private float xMin, xMax, yMin, yMax;

    static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_PREFERRED;

    public Mandelbrot() {
	
    }

    public Mandelbrot(float maxIter, float maxSquareModulus) {
	this.maxIter = maxIter;
	this.maxSquareModulus = maxSquareModulus;
    }

    public Mandelbrot(float[] params) {
	setAll(params);
    }

    public float getMaxIter() {
	return maxIter;
    }

    // set the region to be considered from points
    public void setRegion(float xMin, float xMax, float yMin, float yMax) {
	this.xMin = xMin;
	this.xMax = xMax;
	this.yMin = yMin;
	this.yMax = yMax;
    }

    // set the region to be considered from an array of coordinates
    public void setRegion(float[] coords) {
	this.xMin = coords[0];
	this.xMax = coords[1];
	this.yMin = coords[2];
	this.yMax = coords[3];
    }

    public void setIterAndModulus(float maxIter, float maxSquaredModulus) {
	this.maxIter = maxIter;
	this.maxSquareModulus = maxSquareModulus;
    }

    // set all parameters; the first four values in params are
    // interpreted as coordinates of the region, while params[4] and
    // params[5] are the maximum iterations and maximum squared
    // modulus, respectively
    public void setAll(float[] params) {
	setRegion(params);
	setIterAndModulus(params[4], params[5]);
    }

    // a baseline implementation of computing escape times for the
    // current region
    // esc is a 2d array to record the escape times,
    // where the first index records rows of the region, and the
    // second index is the column number
    public void escapeTimesBaseline (float[][] esc) {

	float xStep = (xMax - xMin) / esc[0].length;
	float yStep = (yMax - yMin) / esc.length;

    //xMax - xMin is the length of the region we're computing the mandelbrot escape times for
    //yMax - yMin is the height of the region they're computing the escape times for.
    //dividing the region by esc length gives us how far apart each point we're testing is gonna be.

	for (int i = 0; i < esc.length; i++) {
	    for (int j = 0; j < esc[0].length; j++) {
        //we're now gonna work to fill in each entry of esc[i][j] with the escape times of each point (i,j).
		int iter = 0;
		float cx = xMin + j * xStep;
		float cy = yMin + i * yStep;
        //current x and current y of the point that the loop is pointing to.

		float zx = 0;
        //Re(z)
		float zy = 0;
        //Im(z)
		
		while (iter < maxIter && zx * zx + zy * zy < maxSquareModulus) {
		    float z = zx * zx - zy * zy + cx;
            //zx is a and zy is b in a+bi. hence, this reads as
           
		    zy = 2 * zx * zy + cy;
             //Im(z) = 2ab + b
		    zx = z;
             //Re(z) = (a^2) - (b^2) + a
		    iter++;		    
		}

		esc[i][j] = iter;
	    }
	}

    }

    // an optimized implementation of escapeTimesBaseline that uses
    // vector operations
public void escapeTimesOptimized (float[][] esc) {
    
	float xStep = (xMax - xMin) / esc[0].length;
	float yStep = (yMax - yMin) / esc.length;
    //xMax - xMin is the length of the region we're computing the mandelbrot escape times for
    //yMax - yMin is the height of the region they're computing the escape times for.
    //dividing the region by esc length gives us how far apart each point we're testing is gonna be.
	int step = SPECIES.length();
	//int bound1 = SPECIES.loopBound(esc.length);
    int bound2 = SPECIES.loopBound(esc[0].length);  
    //loop bound

    float[] cxa = new float[esc[0].length];
    //array that will store the incremental bits of the cx[] vector to be declared soon
    for(float j=0; j<cxa.length; j=j+1){
        cxa[((int)j)] = j * xStep;
    }
    //adds the corresponding incremental parts to cxa.

    for (float i = 0; i < esc.length; i++) {
    //float cx = xMin + j * xStep;
    
    int j = 0;
    float cy = yMin + i * yStep;
    //cy as defined before
    for (; j < bound2; j += step) {
            //float cx = xMin + j * xStep;                    
            var cx = FloatVector.broadcast(SPECIES, xMin);
            /*if(i == 0 && j == 272){
                    System.out.println("cx is " + cx);
                    System.out.println("cxa[j] is " + cxa[j]);
                }  */
            var cxIncrement = FloatVector.fromArray(SPECIES, cxa, j);
            cx = cx.add(cxIncrement);
            //cx is equal to the incremental part plus the xMin base
             
           /*if(i == 0 && j == 272){
                    System.out.println("cx is " + cx);
                }  
                System.out.println("\n"); */
            
            
            var zx = FloatVector.broadcast(SPECIES, 0);
            var zy = FloatVector.broadcast(SPECIES, 0);
            //set them to 0 as per before
            var iter = FloatVector.broadcast(SPECIES, 0);
            //iter will be what adds into the esc[i][j]

            var bitMask = ((zx.mul(zx)).add(zy.mul(zy))).lt(maxSquareModulus);
            //determines if it's less than the modulus to stop the loop

            int dummyIter = 0;
            while(dummyIter<maxIter && bitMask.anyTrue()) {
            
                bitMask = (((zx.mul(zx)).add(zy.mul(zy))).lt(maxSquareModulus));
                iter = iter.add(1, bitMask);
                //updates the iter to determine if it should be added to an increment of esc[i][j]

                var z = ((zx.mul(zx)).sub(zy.mul(zy))).add(cx);
                //zx is a and zy is b in a+bi. hence, this reads as
                zy = (zx.mul(zy).mul(2)).add(cy);
                //Im(z) = 2ab + b
                zx = z;
                dummyIter++;
                //increments dummyIter
		    }
            iter.intoArray(esc[(int)i], j);
	}
    
    	for (; j < esc[0].length; j++) {
        //we're now gonna work to fill in each entry of esc[i][j] with the escape times of each point (i,j).
        //with every part that steps did not cover above
		int iter = 0;
		float cx = xMin + j * xStep;
        //current x and current y of the point that the loop is pointing to.
        //cy already defined!
		float zx = 0;
        //Re(z)
		float zy = 0;
        //Im(z)
		
		while (iter < maxIter && zx * zx + zy * zy < maxSquareModulus) {
		    float z = zx * zx - zy * zy + cx;
            //zx is a and zy is b in a+bi. hence, this reads as
           
		    zy = 2 * zx * zy + cy;
             //Im(z) = 2ab + b
		    zx = z;
             //Re(z) = (a^2) - (b^2) + a
		    iter++;		    
		}

		esc[(int)i][j] = iter;
        //completes the end part w/old code reuse for once threads are no longer being used.
	    } 
    }
}

}
