package rna1;

import java.util.*;

/**
 * <p>Title: ReadFile</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Luda
 * @version 1.0
 */

//===========================================
// klass petel
//===========================================
class transit_loop
{
  int input_helix; // index vhodyaszhei spirali v massive heliset
  pair input_pair; // para, zakryvayuschaya petlyu snizhu
  Vector output_helix = new Vector(); // indeksy vyhodyaszhih spiralei
  Vector output_pairs = new Vector(); // vektor par-koncov vyhodyaszhih spiralei
  double energy;
  heliset heliset1; // spirali

  //=========================================
  public transit_loop(pair in_pair, int in_helix, Vector out_pairs, Vector out_helix, heliset input_heliset)
  {
    input_pair = new pair(in_pair.l, in_pair.r);
    input_helix = in_helix;
    for (int i = 0; i < out_pairs.size(); i++)
    {
      output_pairs.add(out_pairs.get(i));
      output_helix.add(out_helix.get(i));
    }
    heliset1 = input_heliset;
    energy_calculation();
  }

//=========================================
  void energy_calculation()
  {
   // helix in_helix = (helix) heliset1.get(input_helix); // vhodyaschaya spiral'
    int loop_length = 0;
    ////////////////////////////////
    // kogda shpilka
    if (output_helix.size() == 0)
    {
     //energy = 0;
     loop_length = input_pair.r - input_pair.l-1; // dlina petli
     hairpin_energy(loop_length);//*/
    }
    ///////////////////////////////////////
    // kogda vnutrennyaya petlya ili vypyachivanie
    if (output_helix.size() == 1)
    {
      //energy = 0;
      Integer out_index = (Integer) output_helix.get(0); // index vyhodyaschei spirali
      pair out_pairs = (pair) output_pairs.get(0); // vyhodyaschaya para
  //    if (((out_helix.ls1-in_helix.ls2-1) > 2)&&((in_helix.rs1-out_helix.rs2-1)>2)) // kogda petlya, a ne vypyachivanie
     int l_loop = out_pairs.l-input_pair.l-1; // dlina petli na levom pleche
     int r_loop = input_pair.r-out_pairs.r-1; // dlina petli na pravom pleche
     if (l_loop < 0) l_loop = 0;
     if (r_loop < 0) r_loop = 0;
     loop_length = l_loop + r_loop; // dlina vnutrennei petli
     // proveryaem, chto eti spirali ne obrazuyut skol. spiral'
     int slid_index = heliset1.Compatible[input_helix][out_index.intValue()];
     if (slid_index > 0)
     {
       // esli obrazuet, to dobavlyaem dlinu skol. uchactka, kotoryi raven peresecheniyu plechei
       slid_helix sh = (slid_helix) heliset1.slid_heliset[slid_index];
       loop_length += sh.cross;
     }
     internal_energy(loop_length); //*/
   }
    //////////////////////////////////////////////
    // kogda multipetlya
    if (output_helix.size() > 1)
    {
      //energy = 0;
      loop_length = 0;
      this.sort(); // uporaydochivaem output_helix po vozrastaniyu
      int in_index = input_helix; // nomer vhodyaschei spirali
      int in_B = input_pair.l; // koordinaty vhodnoi spirali
      int in_C = input_pair.r;
      int output_count = output_helix.size(); // kolichestvo vyhodyaschih spiralei
      for (int i = 0; i < output_count; i++)
      {
        Integer out_index = (Integer) output_helix.get(i); // index vyhodyaschei spirali
        pair out_pair = (pair) output_pairs.get(i); // vyhodyaschaya spiral'
        int out_A = out_pair.l;
        int out_D = out_pair.r;
        if (out_A > in_B)
        {
          loop_length += (out_A - in_B - 1);
        }
        // proveryaem, chto spirali in_index i out_index ne obrazuyut skol. spiral'
        int slid_index = heliset1.Compatible[in_index][out_index.intValue()];
        if (slid_index > 0)
        {
          // esli obrazuet, to dobavlyaem dlinu skol. uchactka, kotoryi raven peresecheniyu plechei
          slid_helix sh = (slid_helix) heliset1.slid_heliset[slid_index];
          loop_length += sh.cross;
        }
        // sdvigaem na uzhe obrabotannyi kusok petli
        in_B = out_D;
        in_index = out_index.intValue();
      }
      if ( (in_C - in_B - 1) >= 0)
        loop_length += (in_C - in_B - 1);
      int slid_index = heliset1.Compatible[in_index][input_helix];
      if (slid_index > 0)
      {
        // esli obrazuet, to dobavlyaem dlinu skol. uchactka, kotoryi raven peresecheniyu plechei
        slid_helix sh = (slid_helix) heliset1.slid_heliset[slid_index];
        loop_length += sh.cross;
      }
      loop_length += output_count; // plus kolichestvo spiralei v petle
      internal_energy(loop_length);// schitaem energiyu vnutrenei petli//*/
    }
  }

  // vychislyaem energiyu spilki
  void hairpin_energy(int loop_length)
  {
    if (loop_length == 3)
      energy = 5.7*rna1.energy_loops_coef[0];
    if ((loop_length == 4) || (loop_length == 5)|| (loop_length == 8))
      energy = 5.6*rna1.energy_loops_coef[1];
    if (loop_length == 6)
      energy = 5.4*rna1.energy_loops_coef[2];
    if (loop_length == 7)
      energy = 5.9*rna1.energy_loops_coef[3];
    if (loop_length == 9)
      energy = 6.4*rna1.energy_loops_coef[4];
    if (loop_length == 10)
      energy = 6.5*rna1.energy_loops_coef[5];
    if (loop_length == 11)
      energy = 6.6*rna1.energy_loops_coef[6];
    if (loop_length == 12)
      energy = 6.7*rna1.energy_loops_coef[7];
    if (loop_length == 13)
      energy = 6.8*rna1.energy_loops_coef[8];
    if ((loop_length == 14) || (loop_length == 15))
      energy = 6.9*rna1.energy_loops_coef[9];
    if (loop_length == 16)
      energy = 7.0*rna1.energy_loops_coef[10];
    if ((loop_length == 17) || (loop_length == 18))
      energy = 7.1*rna1.energy_loops_coef[11];
    if ((loop_length == 19) || (loop_length == 20))
      energy = 7.2*rna1.energy_loops_coef[12];
    if ((loop_length == 21) || (loop_length == 22))
      energy = 7.3*rna1.energy_loops_coef[13];
    if ((loop_length == 23) || (loop_length == 24))
      energy = 7.4*rna1.energy_loops_coef[14];
    if ((loop_length == 25) || (loop_length == 26)|| (loop_length == 27))
      energy = 7.5*rna1.energy_loops_coef[15];
    if ((loop_length == 28) || (loop_length == 29))
      energy = 7.6*rna1.energy_loops_coef[16];
    if (loop_length >= 30) // nado naity formulu dlya petli >30
      energy = 7.7*rna1.energy_loops_coef[17];
    energy = energy*1000; // energiya v kal/mol
  }

  // vychislyaem energiyu vnutrennei petli
  void internal_energy(int loop_length)
  {
    if (loop_length <= 4)
      energy = 1.7*rna1.energy_loops_coef[18];
    if (loop_length == 5)
      energy = 1.8*rna1.energy_loops_coef[19];
    if (loop_length == 6)
      energy = 2.0*rna1.energy_loops_coef[20];
    if (loop_length == 7)
      energy = 2.2*rna1.energy_loops_coef[21];
    if (loop_length == 8)
      energy = 2.3*rna1.energy_loops_coef[22];
    if (loop_length == 9)
      energy = 2.4*rna1.energy_loops_coef[23];
    if (loop_length == 10)
      energy = 2.5*rna1.energy_loops_coef[24];
    if (loop_length == 11)
      energy = 2.6*rna1.energy_loops_coef[25];
    if (loop_length == 12)
      energy = 2.7*rna1.energy_loops_coef[26];
    if (loop_length == 13)
      energy = 2.8*rna1.energy_loops_coef[27];
    if (loop_length == 14)
      energy = 2.9*rna1.energy_loops_coef[28];
    if ((loop_length == 15) || (loop_length == 16))
      energy = 3.0*rna1.energy_loops_coef[29];
    if ((loop_length == 17) || (loop_length == 18))
      energy = 3.1*rna1.energy_loops_coef[30];
    if (loop_length == 19)
      energy = 3.2*rna1.energy_loops_coef[31];
    if ((loop_length == 20) || (loop_length == 21))
      energy = 3.3*rna1.energy_loops_coef[32];
    if ((loop_length == 22) || (loop_length == 23) || (loop_length == 24))
      energy = 3.4*rna1.energy_loops_coef[33];
    if ((loop_length == 25) || (loop_length == 26))
      energy = 3.5*rna1.energy_loops_coef[34];
    if ((loop_length == 27) || (loop_length == 28) || (loop_length == 29))
      energy = 3.6*rna1.energy_loops_coef[35];
    if (loop_length >= 30) // nado naity formulu dlya petli >30
      energy = 3.7*rna1.energy_loops_coef[36];
    energy = energy*1000; // energiya v kal/mol
  }

  // uporyadochivaem vector output_helix
  void sort()
  {
    for (int i = 0; i < output_helix.size(); i++)
    {
      Integer sh = (Integer) output_helix.get(i);
      for (int j = i + 1; j < output_helix.size(); j++)
      {
        Integer sh1 = (Integer) output_helix.get(j);
        pair pair1 = (pair) output_pairs.get(j);
        if (sh.intValue() > sh1.intValue())
        {
          output_helix.remove(j);
          output_pairs.remove(j);
          output_helix.add(i, sh1);
          output_pairs.add(i, pair1);
        }
      }
    }
    int t = -1;
    for (int i = 0; i < output_helix.size(); i++)
    {
      Integer sh = (Integer) output_helix.get(i);
      int ls = sh.intValue();
      for (int j = i + 1; j < output_helix.size(); j++)
      {
        Integer sh1 = (Integer) output_helix.get(j);
        if (sh1.intValue() < ls)
        {
          ls = sh1.intValue();
          t = j;
        }
      }
      if (t != -1)
      {
        Integer sh1 = (Integer) output_helix.get(t);
        pair pair1 = (pair) output_pairs.get(t);
        output_helix.remove(t);
        output_pairs.remove(t);
        output_helix.add(i, sh1);
        output_pairs.add(i, pair1);
        t = -1;
      }
      else continue;
    }
  }

}

//================================================
// klass sparennyh ocnovanii (pary)
//================================================
class pair
{
  int l;
  int r;

  public pair (int x1, int y1)
  {
    l = x1;
    r = y1;
  }
}


//==================================================
// klass perehodnogo sostoyaniya
//==================================================
public class transit_foldings /*extends rna1.foldings*/
{
  private int[] down_loop; // index petli v folding, iz kotoroi spiral vyhodit
  private int[] up_loop; // index petli v folding, v kotoruyu spiral vhodit
  private heliset heliset1;
  int heliset_size;

  Vector folding = new Vector(); // structura RNK (derevo), vektor sostoit iz ob'ektov klassa loop
  //double folding_energy = 0; // energiya etoi structury
  double loops_energy = 0; // energiya petel' v etoi sturcture
  boolean exists[];
  pair[] pairsAD; // vektor par AD
  pair[] pairsBC; // vektor par BC

  public transit_foldings(pair[] pairsAD1, pair[] pairsBC1, boolean exists1[], heliset in_heliset)
  {
    heliset1 = in_heliset;
    heliset_size = heliset1.heliset_size;
    exists = new boolean[heliset_size];
    pairsAD = new pair[heliset_size];
    pairsBC = new pair[heliset_size];
    for (int i = 0; i < heliset_size; i++)
    {
      exists[i] = exists1[i];
      pairsAD[i] = pairsAD1[i];
      pairsBC[i] = pairsBC1[i];
    }
    down_loop = new int[heliset_size]; // index petli v folding, iz kotoroi spiral vyhodit
    up_loop = new int[heliset_size]; // index petli v folding, v kotoruyu spiral vhodit

    loops_add();
  }

  //====================================
  // kopiruem ob'ekt

  public transit_foldings(transit_foldings f)
  {
    heliset1 = f.heliset1;
    heliset_size = heliset1.heliset_size;
    folding = (Vector) f.folding.clone();
    down_loop = new int[heliset_size]; // index petli v folding, iz kotoroi spiral vyhodit
    up_loop = new int[heliset_size]; // index petli v folding, v kotoruyu spiral vhodit
    for (int i = 0; i < heliset_size; i++)
    {
      down_loop[i] = f.down_loop[i];
      up_loop[i] = f.up_loop[i];
    }
    exists = new boolean[heliset_size];
    for (int i = 0; i < heliset_size; i++)
    {
      exists[i] = f.exists[i];
    }
    loops_energy = f.loops_energy;
  }

  //=======================================
  // schitaem energiyu structury folding
  void energy_calculation()
  {
    loops_energy = 0;
 //   double helises_energy = 0;
    // schitaem energiyu petel'
    for (int i = 0; i < folding.size(); i++)
    {
      transit_loop loop2 = (transit_loop) folding.get(i);
      loops_energy += loop2.energy;
      // schitaem energiyu spiralei
  //    helix helix1 = (helix) heliset1.get(loop2.input_helix); // spiral
  //    helises_energy += helix1.energy; // dobavlyaem energiyu spirali
    }
    // vychitaem energiyu peresecheniya dlya skol. spiralei
/*    Vector exists_helises = index_exists_true(); // nomera suschestvuyuschih spiralei
    for (int i = 0; i < exists_helises.size(); i++)
    {
      int I = ( (Integer) exists_helises.get(i)).intValue();
      for (int j = i; j < exists_helises.size(); j++)
      { // esli eta spiral' suchestvuet
        int J = ( (Integer) exists_helises.get(j)).intValue();
        if (heliset1.Compatible[I][J] > 0) // esli skol. spiral'
        //I>J is added by A.Favorov
        {
          slid_helix sh = (slid_helix) heliset1.slid_heliset.get(heliset1.Compatible[I][J]);
          // umenishaem energiyu spiralei ne energiyu peresecheniya
          helises_energy -= sh.energy_cross;
        }
      }
    }*/
    // schitaem energiyu struktuy
   // folding_energy = loops_energy + helises_energy;
  }

  void loops_add()
  {
    folding.removeAllElements();
    for (int i = 0; i < heliset_size; i++)
    {
      down_loop[i] = -1;
      up_loop[i] = -1;
    }
    for (int k = 0; k < heliset_size; k++)
    {
      if (exists[k])
      {
        Vector upper = upper_helix(k);
        Vector upper_pairs = new Vector();
        for (int i = 0; i < upper.size(); i++)
        {
          Integer t = (Integer) upper.get(i);
          upper_pairs.add(i, pairsAD[t.intValue()]);
        }
        pair k_pair = pairsBC[k];
        folding.add(new transit_loop(k_pair, k, upper_pairs, upper, heliset1)); // dobavlyaem petlyu v folding
        int loop_index = folding.size() - 1;
        // zapolnyaem massivy up_loop i down_loop
        up_loop[k] = loop_index;
        for (int i = 0; i < upper.size(); i++)
        {
          Integer k1 = (Integer) upper.get(i);
          down_loop[k1.intValue()] = loop_index;
        }
      }
    }
    energy_calculation();
  }

  //==================================
  // iszhem spirali vyshe k-oi
  Vector upper_helix(int k)
  {
    Vector upper = new Vector(); // shpilki vyshe rassmatrivaemoi
    pair sh_k = pairsBC[k];
    int l = sh_k.l;
    int r = sh_k.r;

    for (int i = heliset_size-1; i >=0; i--)
    {
      if (exists[i]) // esli sovmestimy i mogut suszhestvovat'
      {
        pair sh_i = pairsAD[i];
        if ( (sh_i.l > l) && (sh_i.r < r))
        {
          upper.add(new Integer(i));
          r = sh_i.r;
        }
      }
    }
    return upper;
  }


  //===================================
  // dobavit' t-uyu spiral' v folding
  public void add_helix(int t)
  {
    exists[t] = true;
    loops_add();
  }

//===================================
  // udalit' t-uyu spiral' iz folding
  public void delete_helix(int t)
  {
    exists[t] = false;
    loops_add();
  }

  //===================================
  // kak izmenitsya energiya, esli dobavit' t-uyu spiral' v folding
  public double if_add_helix(int t)
  {
    this.add_helix(t);
    double energy = this.loops_energy;
    this.delete_helix(t); // udalyaem tu spiral, kotoruyu dobavili
    return energy;
  }

  //=========================================
  // zapis' struktury v skobochno-tochechnom vide
  String bracket_dots(int sequence_length)
  {
    char[] output_folding = new char[sequence_length]; // massiv skobochek i tochechek
    for (int i = 0; i < sequence_length; i++) // vezde propisyvaem tochki
    {
      output_folding[i] = '.';
    }
    Vector slid_index = new Vector(); // nomera elementov so skol'zyaschei petlei
    Vector exists_helises = index_exists_true(); // nomera suschestvuyuschih spiralei
    for (int i = 0; i < exists_helises.size(); i++)
    {
      boolean correct_helix = true; // vhodit ili net spiral' v structuru so skol'zyaschei petlei
      int I = ( (Integer) exists_helises.get(i)).intValue();
      for (int j = 0; j < exists_helises.size(); j++)
      { // esli eta spiral' suchestvuet
        int J = ( (Integer) exists_helises.get(j)).intValue();
        if (heliset1.Compatible[I][J] > 0) // pravil'naya spiral' ili net
        {// esli net, to zapominaem nomer structury
          slid_index.add(new Integer(heliset1.Compatible[I][J]));
          correct_helix = false;
        }
      }
      if (correct_helix) // esli pravil'naya spiral', to prosto prostavlyaem skobki
      {
        helix helix1 = (helix) heliset1.get(I); // vytaskivaem sootvetstvuyuschuyu spiral'
        for (int j = helix1.A; j <= helix1.B; j++)
        { // levoe plecho
          output_folding[j] = '(';
        }
        for (int j = helix1.C; j <= helix1.D; j++)
        { // pravoe plecho
          output_folding[j] = ')';
        }
      }
    }
    StringBuffer slid_helises = new StringBuffer();
    // esli est' skol'zyaschie petli
    if (!slid_index.isEmpty())
    {
      // ubiraem odinakovye
      for (int i = 0; i < slid_index.size(); i++)
      {
        Integer I = (Integer) slid_index.get(i);
        for (int j = i + 1; j < slid_index.size(); j++)
        {
          Integer J = (Integer) slid_index.get(j);
          if (I.intValue() == J.intValue())
          {
            slid_index.remove(j--);
          }
        }
      }
      // zapolnyaem pro skol'zyaschie petli
      for (int i = 0; i < slid_index.size(); i++)
      {
        int index = ( (Integer) slid_index.get(i)).intValue();
        slid_helix helix1 = (slid_helix) heliset1.slid_heliset[index];
        slid_helises.append(helix1.toString()); // zapisali peresekayuschiesya pozicii
        for (int j = helix1.A1; j <= helix1.B1; j++)
        { // levoe plecho
          output_folding[j] = '(';
        }
        for (int j = helix1.C1; j <= helix1.D1; j++)
        { // pravoe plecho
          output_folding[j] = ')';
        }
        for (int j = helix1.A2; j <= helix1.B2; j++)
        { // levoe plecho
          output_folding[j] = '(';
        }
        for (int j = helix1.C2; j <= helix1.D2; j++)
        { // pravoe plecho
          output_folding[j] = ')';
        }
      }
    }
    StringBuffer sb = new StringBuffer(sequence_length);
    for (int i = 0; i< sequence_length; i++)
    {
      sb.append(output_folding[i]);
    }
    if (!slid_index.isEmpty()) sb.append("\t" + slid_helises.toString()); // dobavili peresekayuschiesya pozicii posle skobok-tochek cherez tab
    return sb.toString();
  }

  //=================================
  // vydaet vektor s nomerami elementov v exists, kotorye ravny true, t.e. sootvetstvuyuschie spirali suschestvuyut
  Vector index_exists_true()
  {
    Vector true_elem = new Vector(); // vybiraem elementy ravnye true
    for (int j = 0; j < exists.length; j++)
    {
      if (exists[j] == true)
      {
        true_elem.add(new Integer(j));
      }
    }
    return true_elem;
  }

/*//=====================================
  public String toString(int length)
  {
    StringBuffer sb = new StringBuffer();
    sb.append("structure energy = " + folding_energy +"\n" +
              "loops energy = " + loops_energy + "\n" +
               bracket_dots(length) + "\n");
    return sb.toString();
  }

  public String print_helises(boolean only_index)
  {
    StringBuffer sb = new StringBuffer();
    Vector true_elem = index_exists_true();
    // pechataem vybrannye elementy, t.e. suschestvuyuschie spirali
    for (int j = 0; j < true_elem.size(); j++)
    {
      Integer J = (Integer) true_elem.get(j);
      sb.append(J + "\n");
      if (!only_index)
        {
          helix helix1 = (helix) heliset1.get(J.intValue());
          sb.append(helix1.toString(false));
        }
    }
    return sb.toString();
  }*/

}



