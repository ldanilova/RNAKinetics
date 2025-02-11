package rna1;

//import java.util.*;

/**
 * <p>Title: ReadFile</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Luda
 * @version 1.0
 */

public class ArrayRandom {

  int draw(double weights[], double k)
  {
    // modeliruem sluchainuyu velichinu r, ravnomerno respredelennuyu na [0,k]

/*   Random Random1 = new Random(150); //*/
    double g = rna1.Random1.nextDouble();
 //  double g = Math.random();

    double r = k*g;
    // nado naiti takoe i, chto sum kj (j=1:i-1) <= r < sum kj (j=1:i)
    int i_per=-1; // to samoe i
    double bound_up = weights[0]; // sum kj (j=1:i-1)
    for (int i = 1; i <= weights.length; i++)
    {
      if (r<bound_up)
      {
        i_per=i-1;
        break;
      }
      bound_up = bound_up+weights[i];
    }
    return i_per;
  }

  int draw(double weights[])
  {
    double k=0;
    for (int i = 0; i < weights.length; i++) k = k + weights[i];
    // modeliruem sluchainuyu velichinu r, ravnomerno respredelennuyu na [0,k]


  /*  Random Random1 = new Random(150);
    double g = Random1.nextDouble(); //*/


    double g = Math.random();

    double r = k*g;
    // nado naiti takoe i, chto sum kj (j=1:i-1) <= r < sum kj (j=1:i)
    int i_per=-1; // to samoe i
    double bound_up = weights[0]; // sum kj (j=1:i-1)
    for (int i = 1; i <= weights.length; i++)
    {
      if (r<bound_up)
      {
        i_per=i-1;
        break;
      }
      bound_up = bound_up+weights[i];
    }
    return i_per;
  }
}
