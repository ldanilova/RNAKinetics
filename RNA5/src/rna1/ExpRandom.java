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

public class ExpRandom
{
  double draw(double k)// vydaet sluchainuyu velichinu
    {
    /* Random Random1 = new Random(150);//*/
      double g = rna1.Random1.nextDouble();
     // double g = Math.random();

      double tau = - (1.0/k)*Math.log(g);
      return tau;
  }
}
