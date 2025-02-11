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

//$Id: helix.java,v 1.11 2003/10/29 01:35:55 favorov Exp $

class helix {
  int A, B; // nachalo i konec levogo plecha
  int C, D; // nachalo i konec pravogo plecha
  String Left;
  String Right; // reversed rigth strand
  double energy;// energiya shpilki v kal/mol'
  double k_break; // konstanta skorosti razrusheniya
  double[] graph; // grafik
  double grown_time; // moment vremeni, kogda eta spiral' vyrosla

   // konstanty dlyz podscheta konstanty razrusheniya
  final double k0 = rna1.k0; // konstanta skorosti obrazovaniya ocherednoi pary (c^-1)
  final double R = 1.987; //gazovaya postoyannaya 8.314 Dj/K*mol' (1.987 kal/K*mol')
  final double T = 310; // temperature v K
  final double RT = R*T;

// konstructor dlya spirali bez skol'zyaschei petli
   public helix(int A_par, int B_par, int C_par, int D_par, String seq)
   {
     A = A_par;
     B = B_par;
     C = C_par;
     D = D_par;
     Left = seq.substring(A, B + 1); // levoe plecho
     Right = reverse(seq.substring(C, D + 1)); // perevernutoe pravoe
     energy = energy_calculation(Left, Right);
     k_break = ( (Right.length()-1) * k0 * Math.exp(energy / (RT)));
     graph = new double[rna1.graph_step];
     grown_time = 0;
   };

   //==================================
  // perevorachivaet stroku s
  String reverse(String s)
  {
    String s1 = "";
    for (int i = s.length() - 1; i >= 0; i--)
    {
      s1 = s1.concat(s.substring(i, i + 1));
    }
    return s1;
  }

   //=======================================
   // schitaem energiyu shpilki, kogda energiya zavisit ot predyduschei pary
   double energy_calculation(String Left, String Right)
   {
     double en=0;
     for (int i=0; i < Left.length()-1; i++)
     {
       if ( (Left.charAt(i) == 'a') && (Right.charAt(i) == 'u'))
       {
         if ( (Left.charAt(i + 1) == 'a') && (Right.charAt(i + 1) == 'u'))
           en = en - 0.9*rna1.energy_helix_coef[0];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'a'))
           en = en - 1.1*rna1.energy_helix_coef[1];
         if ( (Left.charAt(i + 1) == 'c') && (Right.charAt(i + 1) == 'g'))
           en = en - 2.2*rna1.energy_helix_coef[2];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'c'))
           en = en - 2.1*rna1.energy_helix_coef[3];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'g'))
           en = en - 1.4*rna1.energy_helix_coef[4];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'u'))
           en = en - 0.6*rna1.energy_helix_coef[5];
           continue;
       }
       if ( (Left.charAt(i) == 'u') && (Right.charAt(i) == 'a'))
       {
         if ( (Left.charAt(i + 1) == 'a') && (Right.charAt(i + 1) == 'u'))
           en = en - 1.3*rna1.energy_helix_coef[6];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'a'))
           en = en - 0.9*rna1.energy_helix_coef[7];
         if ( (Left.charAt(i + 1) == 'c') && (Right.charAt(i + 1) == 'g'))
           en = en - 2.4*rna1.energy_helix_coef[8];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'c'))
           en = en - 2.1*rna1.energy_helix_coef[9];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'g'))
           en = en - 1.3*rna1.energy_helix_coef[10];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'u'))
           en = en - 1.0*rna1.energy_helix_coef[11];
           continue;
       }
       if ( (Left.charAt(i) == 'g') && (Right.charAt(i) == 'c'))
       {
         if ( (Left.charAt(i + 1) == 'a') && (Right.charAt(i + 1) == 'u'))
           en = en - 2.4*rna1.energy_helix_coef[12];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'a'))
           en = en - 2.2*rna1.energy_helix_coef[13];
         if ( (Left.charAt(i + 1) == 'c') && (Right.charAt(i + 1) == 'g'))
           en = en - 3.4*rna1.energy_helix_coef[14];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'c'))
           en = en - 3.3*rna1.energy_helix_coef[15];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'g'))
           en = en - 2.5*rna1.energy_helix_coef[16];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'u'))
           en = en - 1.5*rna1.energy_helix_coef[17];
           continue;
       }
       if ( (Left.charAt(i) == 'c') && (Right.charAt(i) == 'g'))
       {
         if ( (Left.charAt(i + 1) == 'a') && (Right.charAt(i + 1) == 'u'))
           en = en - 2.1*rna1.energy_helix_coef[18];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'a'))
           en = en - 2.1*rna1.energy_helix_coef[19];
         if ( (Left.charAt(i + 1) == 'c') && (Right.charAt(i + 1) == 'g'))
           en = en - 3.3*rna1.energy_helix_coef[20];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'c'))
           en = en - 2.4*rna1.energy_helix_coef[21];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'g'))
           en = en - 2.1*rna1.energy_helix_coef[22];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'u'))
           en = en - 1.4*rna1.energy_helix_coef[23];
           continue;
       }
       if ( (Left.charAt(i) == 'g') && (Right.charAt(i) == 'u'))
       {
         if ( (Left.charAt(i + 1) == 'a') && (Right.charAt(i + 1) == 'u'))
           en = en - 1.3*rna1.energy_helix_coef[24];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'a'))
           en = en - 1.4*rna1.energy_helix_coef[25];
         if ( (Left.charAt(i + 1) == 'c') && (Right.charAt(i + 1) == 'g'))
           en = en - 2.5*rna1.energy_helix_coef[26];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'c'))
           en = en - 2.1*rna1.energy_helix_coef[27];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'g'))
           en = en + 1.3*rna1.energy_helix_coef[28];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'u'))
           en = en - 0.5*rna1.energy_helix_coef[29];
           continue;
       }
       if ( (Left.charAt(i) == 'u') && (Right.charAt(i) == 'g'))
       {
         if ( (Left.charAt(i + 1) == 'a') && (Right.charAt(i + 1) == 'u'))
           en = en - 1.0*rna1.energy_helix_coef[30];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'a'))
           en = en - 0.6*rna1.energy_helix_coef[31];
         if ( (Left.charAt(i + 1) == 'c') && (Right.charAt(i + 1) == 'g'))
           en = en - 1.5*rna1.energy_helix_coef[32];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'c'))
           en = en - 1.4*rna1.energy_helix_coef[33];
         if ( (Left.charAt(i + 1) == 'u') && (Right.charAt(i + 1) == 'g'))
           en = en - 0.5*rna1.energy_helix_coef[34];
         if ( (Left.charAt(i + 1) == 'g') && (Right.charAt(i + 1) == 'u'))
           en = en + 0.3*rna1.energy_helix_coef[35];
           continue;
       }
     }
     return (en*1000); // energiya v kal/mol
    }

    //===============================
  public String toString()
     {
       StringBuffer sb = new StringBuffer();
       sb.append("" + A + "-" + B + ":" + D + "-" + C + "\n" +
               Left + "\n" + Right + "\n" +
               "helix energy = "+ energy + "\n" +
               "helix length = " + Left.length() + "\n" +
               "k_break = " + k_break + "\n"/* +
               "grown time = " + grown_time + "\n"*/);
     if (rna1.print_graph == true)
     {
       sb.append("graph 2: ");
       for (int i = 0; i < rna1.graph_step; i++)
       {
         // delim na kolichestvo prohodov
         graph[i] = (graph[i] / (rna1.M));
         // delim na shag, chtoby poluchilas' veroyatnost'
         graph[i] = (graph[i] / (rna1.Teta/rna1.graph_step));
         sb.append(graph[i] + ", ");
       }
       sb.append("\n");
     }
       return sb.toString();
     };


  public String toString(boolean graphic)
  {
    StringBuffer sb = new StringBuffer();
    sb.append("" + A + "-" + B + ":" + D + "-" + C + "\n" +
              Left + "\n" + Right + "\n" +
              "helix energy = "+ energy + "\n" +
              "helix length = " + Left.length() + "\n" +
              "k_break = " + k_break + "\n");
    if (graphic == true)
    {
      sb.append("graph 2: ");
      for (int i = 0; i < rna1.graph_step; i++) {
        graph[i] = (graph[i] / (rna1.M));
        sb.append(graph[i] + ", ");
      }
    }
    sb.append("\n");
    return sb.toString();
  };

  /*  //============================
  // schitaem energiyu shpilki
  double energy_calculation() {
    double en = 0;
    for (int i = 1; i <= this.Left.length() - 1; i++)
    {
      char lch = Left.charAt(i);
      char rch = Right.charAt(i);
      if ( ( (lch == 'c') && (rch == 'g')) || ( (lch == 'g') && (rch == 'c')))
        en = en - 3;
      if ( ( (lch == 'a') && (rch == 'u')) || ( (lch == 'u') && (rch == 'a')))
        en = en - 2;
      if ( ( (lch == 'g') && (rch == 'u')) || ( (lch == 'u') && (rch == 'g')))
        en = en - 1;
    }
    return en;
  }*/

}
//=====================================================
// klass par spiralei so skol'zyaschei petlei
class slid_helix
{
  // sleduyuschie parametry otnosyatsya k spiralyam so skol'zyaschei petlei
  int slip_parts[]; // massiv, kotoryi govorit kakie spirali vhodyat v nee
  // energii kuska ot pervoi i vtoroi spiralei, ih peresecheniya i sk. petli
  double energy1, energy2, energy_cross;
  int cross = 0;//dlina peresecheniya
  int cross_pos[]; // nomera peresecheniya
  int pos_h1[]; // nomera al'ternativnyh pozicii iz spirali 1
  int pos_h2[]; // nomera al'ternativnyh pozicii iz spirali 2
  int A1, B1, C1, D1, A2, B2, C2, D2; // pravye i levye plechi spiralei bez peresecheniya
  boolean left_slid = false; // skol'zyaschaya petlya na levom pleche
  boolean bottom_slid = false; // spirali peresekayutsya neodnoimennymi plechami

// konstructor dlya spirali so skol'zyaschei petlei, kotoraya poluchikas' iz spiralei h1 i h2
   public slid_helix(int h1, int h2, helix helix1, helix helix2, boolean bottom)
   {
     slip_parts = new int[2];
     slip_parts[0] = h1;
     slip_parts[1] = h2;
     int n1 = helix1.Left.length(); // dlina suschestvuyuschei spirali
     int n2 = helix2.Left.length(); // dlina spirali, kotoraya hochet obrazovat'sya
     int A_1 = helix1.A; // nachalo levogo plecha
     int A_2 = helix2.A;
     int B_1 = helix1.B; // konec levogo plecha
     int B_2 = helix2.B;
     int C_1 = helix1.C; // nachalo pravogo plecha
     int C_2 = helix2.C;
     int D_1 = helix1.D; // konec pravogo plecha
     int D_2 = helix2.D;
     if (!bottom) // peresekayutsya odnoimennye plechi
     {
       if (A_1 < A_2) // helix1 lezhit nizhe
       {
         cross = 0;
         if ((D_2 >= C_1)&&( A_2 > B_1 )) // peresechenie pravyh plechei, petlya na levom, levye ne peresekayutsya
         {
           left_slid = true;
           cross = D_2 - C_1 + 1; // peresechenie pravyh plechei
           cross_pos = new int[cross];
           pos_h1 = new int[cross];
           pos_h2 = new int[cross];
           for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
           {
             cross_pos[i] = C_1 + i;
             pos_h1[i] = B_1 - i;
             pos_h2[i] = A_2 + cross - 1 - i;
           }

         }
         if (D_2 < C_1 && A_2 <= B_1) // peresecheniy levyh plechei, pelya na pravom pleche, pravye ne peresekayutsya
         {
           cross = B_1 - A_2 + 1; // peresechenie levyh plechei
           cross_pos = new int[cross];
           pos_h1 = new int[cross];
           pos_h2 = new int[cross];
           for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
           {
             cross_pos[i] = A_2 + i;
             pos_h1[i] = C_1 + cross - 1 - i;
             pos_h2[i] = D_2 - i;
           }

         }
         if (A_2 <= B_1 && D_2 >= C_1) // peresekayutsya oba plecha
         {
           int cross_l = B_1 - A_2 + 1; // peresechenie levyh plechei
           int cross_r = D_2 - C_1 + 1; // peresechenie pravyh plechei
           if (cross_l > cross_r) // esli bol'she peresechenie na levyh plechah
             {
               cross = cross_l;
               cross_pos = new int[cross];
               pos_h1 = new int[cross];
               pos_h2 = new int[cross];
               for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
               {
                 cross_pos[i] = C_1 + i;
                 pos_h1[i] = B_1 - i;
                 pos_h2[i] = A_2 + cross - 1 - i;
               }

             }
           else // esli bol'she peresechenie na pravyh plechah
           {
             cross = cross_r;
             left_slid = true;
             cross_pos = new int[cross];
             pos_h1 = new int[cross];
             pos_h2 = new int[cross];
             for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
             {
               cross_pos[i] = A_2 + i;
               pos_h1[i] = C_1 + cross - 1 - i;
               pos_h2[i] = D_2 - i;
             }
           }
         }
         A1 = A_1;
         B1 = B_1 - cross;
         C1 = C_1 + cross;
         D1 = D_1;
         A2 = A_2 + cross;
         B2 = B_2;
         C2 = C_2;
         D2 = D_2 - cross;
         // energiya helix1, bez peresecheniya

         energy1 = helix1.energy_calculation(helix1.Left.substring(0, n1 - cross),
                                             helix1.Right.substring(0, n1 - cross));
         // energiya helix2, bez peresecheniya
         energy2 = helix2.energy_calculation(helix2.Left.substring(cross),
                                             helix2.Right.substring(cross));
         // energiya peresecheniya
         /*    energy_cross = helix2.energy_calculation(helix2.Left.substring(0, cross),
          helix2.Right.substring(0, cross));//*/
          // kak srednee mezhdu raznicei energii
         energy_cross = ( (helix1.energy - energy1) + (helix2.energy - energy2)) /2;
       }
       else //(A_1 < A_2) // helix2 lezhit nizhe
       {
         cross = 0;
         if (D_1 >= C_2 && A_1 > B_2) // peresechenie pravyh, petlya na levom pleche, levye ne peresekayutsya
         {
           left_slid = true;
           cross = D_1 - C_2 + 1; // peresechenie pravyh plechei
           cross_pos = new int[cross];
           pos_h1 = new int[cross];
           pos_h2 = new int[cross];
           for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
           {
             cross_pos[i] = C_2 + i;
             pos_h1[i] = helix2.B - i;
             pos_h2[i] = A_1 + cross - 1 - i;
           }
         }
         if (D_1 < C_2 && A_1 <= B_2) // peresechenie levyh, pelya na pravom pleche, pravye ne peresekayutsya
         {
           cross = B_2 - A_1 + 1; // peresechenie levyh plechei
           cross_pos = new int[cross];
           pos_h1 = new int[cross];
           pos_h2 = new int[cross];
           for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
           {
             cross_pos[i] = A_1 + i;
             pos_h1[i] = C_2 + cross - 1 - i;
             pos_h2[i] = D_1 - i;
           }
         }
         if (D_1 >= C_2 && A_1 <= B_2) // peresekayutsya oba plecha
         {
           int cross_r = D_1 - C_2 + 1; // peresechenie pravyh plechei
           int cross_l = B_2 - A_1 + 1; // peresechenie levyh plechei
           if (cross_l > cross_r) // esli bol'she peresechenie na levyh plechah
           {
             cross = cross_l;
             cross_pos = new int[cross];
             pos_h1 = new int[cross];
             pos_h2 = new int[cross];
             for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
             {
               cross_pos[i] = A_1 + i;
               pos_h1[i] = C_2 + cross - 1 - i;
               pos_h2[i] = D_1 - i;
             }
           }
           else // esli bol'she peresechenie na pravyh plechah
           {
             cross = cross_r;
             left_slid = true;
             cross_pos = new int[cross];
             pos_h1 = new int[cross];
             pos_h2 = new int[cross];
             for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
             {
               cross_pos[i] = C_2 + i;
               pos_h1[i] = helix2.B - i;
               pos_h2[i] = A_1 + cross - 1 - i;
             }
           }
         }
         A1 = A_1 + cross;
         B1 = B_1;
         C1 = C_1;
         D1 = D_1 - cross;
         A2 = A_2;
         B2 = B_2 - cross;
         C2 = C_2 + cross;
         D2 = D_2;
         // energiya helix1, bez peresecheniya
         energy1 = helix1.energy_calculation(helix1.Left.substring(cross),
                                             helix1.Right.substring(cross));
         // energiya helix2, bez peresecheniya
         energy2 = helix2.energy_calculation(helix2.Left.substring(0, n2 - cross),
                                             helix2.Right.substring(0, n2 - cross));
         // energiya peresecheniya
         /*    energy_cross = helix1.energy_calculation(helix1.Left.substring(0, cross),
          helix1.Right.substring(0, cross)); //*/
          // kak srednee mezhdu raznicei energii
         energy_cross = ( (helix1.energy - energy1) + (helix2.energy - energy2)) / 2;
       }
     }
     //================================================================

     else // peresekayutsya raznoimennye plechi (bottom = true)
     {
       cross = 0;
       bottom_slid = true;
       if (C_2 < A_1) // peresekayutsya levoe plecho helix1 i pravoe helix2
       {
         cross = D_2 - A_1 + 1; // peresechenie
         cross_pos = new int[cross];
         pos_h1 = new int[cross];
         pos_h2 = new int[cross];
         for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
         {
           cross_pos[i] = A_1 + i;
           pos_h1[i] = D_1 - i;
           pos_h2[i] = A_2 + cross - 1 - i;
         }
         A1 = A_1 + cross;
         B1 = B_1;
         C1 = C_1;
         D1 = D_1 - cross;
         A2 = A_2 + cross;
         B2 = B_2;
         C2 = C_2;
         D2 = D_2 - cross;

         // energiya helix1, bez peresecheniya
         energy1 = helix1.energy_calculation(helix1.Left.substring(cross),
                                             helix1.Right.substring(cross));
         // energiya helix2, bez peresecheniya
         energy2 = helix2.energy_calculation(helix2.Left.substring(cross),
                                             helix2.Right.substring(cross));
         // energiya peresecheniya
         /*    energy_cross = helix2.energy_calculation(helix2.Left.substring(0, cross),
          helix2.Right.substring(0, cross));//*/
          // kak srednee mezhdu raznicei energii
         energy_cross = ( (helix1.energy - energy1) + (helix2.energy - energy2)) / 2;
       }
       else //(C_2 > A_1) // peresekayutsya pravoe plecho helix1 i levoe helix2
       {

         cross = D_1 - A_2 + 1; // peresechenie
         cross_pos = new int[cross];
         pos_h1 = new int[cross];
         pos_h2 = new int[cross];
         for (int i = 0; i < cross; i++) // zapolnyaem massivy s poziciyami
         {
           cross_pos[i] = A_2 + i;
           pos_h1[i] = A_1 + cross - 1 - i;
           pos_h2[i] = D_2 - i;
         }
         A1 = A_1 + cross;
         B1 = B_1;
         C1 = C_1;
         D1 = D_1 - cross;
         A2 = A_2 + cross;
         B2 = B_2;
         C2 = C_2;
         D2 = D_2 - cross;

         // energiya helix1, bez peresecheniya
         energy1 = helix1.energy_calculation(helix1.Left.substring(cross),
                                             helix1.Right.substring(cross));
         // energiya helix2, bez peresecheniya
         energy2 = helix2.energy_calculation(helix2.Left.substring(cross),
                                             helix2.Right.substring(cross));
         // energiya peresecheniya
         /*    energy_cross = helix1.energy_calculation(helix1.Left.substring(0, cross),
          helix1.Right.substring(0, cross)); //*/
          // kak srednee mezhdu raznicei energii
         energy_cross = ( (helix1.energy - energy1) + (helix2.energy - energy2)) /
             2;
       }

     }

 }


   //===================================
   // zapisyvaem pozicii peresecheniya
   public String toString()
   {
     StringBuffer sb = new StringBuffer();
     for (int i = 0; i < cross; i++)
     {
       sb.append(cross_pos[i] + "=" + pos_h1[i] + ":" + pos_h2[i] + ", ");
     }
     return sb.toString();
   }
}



//=================================================

public class heliset extends java.util.Vector
{
  final int minhelix = rna1.minhelix; // razmer minimalnoi spirali
  //added by A.Favorov
  final int minpartialhelix = rna1.minpartialhelix; // razmer minimalnogo kuska spirali
  // torchashego iz skolziashei spirali
  final int minloop = 3; // razmer minimal'noi petli
  //added by A.Favorov
//  final double min_helix_living_time= ;
  // maksimal'naya energiya spirali
  final double min_energy = rna1.helix_threshold;
  // razreshaet/zapreschaet koncevye GU-pary
  final boolean if_GU = rna1.if_GU;

  int Compatible[][]; // massiv sovmestimosti spiralei
  int bind_pairs[][]; // massiv svyazannyh spiralei
  // nesovmestimy - (-1), sovmestimy - 0,
  // chastichno sovmestimy - nomer spirali so skol'z. petlei v vektore slid_heliset
  int heliset_size; // razmer vectora s normal'nymi spiralyami
  private Vector slid_heliset_process; // vector so skol'zyaschimi spiralyami v processe sborki
  slid_helix[] slid_heliset; // massiv so skol'zyaschimi spiralyami posle obrabotki (budem hranit' spirali s 1, a ne 0)
  //=====================================
  // Konstruktor
  public heliset(String s)
  {
    String seq; // tekushaya posledovatel'nost
    seq = s;
    int sl = s.length();
    int[][] match = new int[sl][sl]; // array of matches for letters
    for (int i = 0; i < sl; i++)
    {
      for (int j = sl - 1; j > i; j--)
      {
        char chi = s.charAt(i), chj = s.charAt(j);
         // proveruaem yavlyaetsya li para (chi,chj) kanonicheskoi ili net
        if (((chi == 'a') && (chj == 'u')) || ((chi == 'u') && (chj == 'a')) ||
            ((chi == 'c') && (chj == 'g')) || ((chi == 'g') && (chj == 'c')))
        {
          match[i][j] = match[j][i] = 2;
        }
        if ( ( (chi == 'u') && (chj == 'g')) || ( (chi == 'g') && (chj == 'u')))
        {
          if (if_GU) match[i][j] = match[j][i] = 2; // razreshaem koncevye GU
          else match[i][j] = match[j][i] = 1;// zapreschaem
          // System.out.print(match[i][j]);
        }
      }
      //System.out.print("\n");
    }
    // pechataem matricu match
   /* for (int i = 0; i < sl; i++)
    {
      System.out.print(s.charAt(i));
      for (int j = 0; j < sl; j++)
      {
        System.out.print(match[i][j]);
      }
      System.out.print("\n");
    }//*/

    // ischem obratnye diagonali v match
    int ls1 = -1, ls2 = -1, rs1 = -1, rs2 = -1;
    for (int i = (2 * minhelix - 1); i < sl; i++)// i - stolbec, j - stroka
    {
      int counter = 0; // schitaet 1 i 2
      for (int j = 0; 2 * j <= i; j++) // idem po obratnoi diagonali do osnovnoi
      {
        if (((match[i - j][j] == 0) && (counter == 0)) ||           // net nachala dlya spirali
            ((match[i - j][j] == 1) && (counter == 0)))  continue;  // ili nachinaetsya s GU
        if ((match[i - j][j] == 2) && (counter == 0))
        {
          ls1 = j; // nachalo levogo plecha
          rs2 = i - j; // konec pravogo plecha
          counter++;
          continue;
        }
        else
        if (((match[i - j][j] == 2)||(match[i - j][j] == 1)) && (counter != 0))
        {
          counter++;
          continue;
        }
        if ( (match[i - j][j] == 0) && (counter != 0)) // t.e. match[i][j]==0
        {
          if (counter >= minhelix) // spiral nabralas' bolshe chem minhelix
          {
            if (match[i - j + 1][j - 1] == 2) // i zakanchivaetsya na GC ili AU
            {
              ls2 = j - 1; // konec levogo plecha
              rs1 = i - j + 1; // nachalo pravogo plecha
              this.add(new helix(ls1, ls2, rs1, rs2, seq));
              ls1 = ls2 = rs1 = rs2 = -1;
            }
            else
            { // esli spiral zakanchivaetsya na GU
              counter--;
              for (int t = 2; t < counter; t++)
              {
                if (counter >= minhelix)
                {
                  if (match[i - j + t][j - t] == 2)
                  {
                    ls2 = j - t; // konec levogo plecha
                    rs1 = i - j + t; // nachalo pravogo plecha
                    this.add(new helix(ls1, ls2, rs1, rs2, seq));
                    ls1 = ls2 = rs1 = rs2 = -1;
                    counter = 0;
                  }
                  if (match[i - j + t][j - t] == 1)
                  {
                    counter--;
                    continue;
                  }
                }
                else  break;
              }
            }
          }
          counter = 0;
        }
      }
    }
    // rassmatrivaem kusok matricy pod obratnoi diagonal'yu

    for (int j = 1; j < sl; j++)
    {
      int i = sl - 1;
      int counter = 0;
      for (int k = 0; k < i - j; k++)
      {
        if (((match[i - k][j + k] == 0)&& (counter == 0)) ||         // net nachala dlya spirali
            ((match[i - k][j + k] == 1)&& (counter == 0)))  continue;// ili nachinaetsya s GU
        if ((match[i - k][j + k] == 2)&& (counter == 0))
        {
          ls1 = j + k; // nachalo levogo plecha
          rs2 = i - k; // konec pravogo plecha
          counter++;
          continue;
        }
        else
        if (((match[i - k][j + k] == 2)|| (match[i - k][j + k] == 1))&& (counter != 0))
        {
          counter++;
          continue;
        }
        if ((match[i - k][j + k] == 0)&& (counter != 0)) // t.e. match[i][j]==0
        {
          if (counter >= minhelix)
          {
            if (match[i - k + 1][j + k - 1] == 2)// i zakanchivaetsya na GC ili AU
            {
              ls2 = j + k - 1; // konec levogo plecha
              rs1 = i - k + 1; // nachalo pravogo plecha
              this.add(new helix(ls1, ls2, rs1, rs2, seq));
              ls1 = ls2 = rs1 = rs2 = -1;
            }
            else // esli spiral zakanchivaetsya na GU
            {
              counter--;
              for (int t = 2; t < counter; t++)
              {
                if (counter >= minhelix)
                {
                  if (match[i - k + t][j + k - t] == 2)
                  {
                    ls2 = j + k - t; // konec levogo plecha
                    rs1 = i - k + t; // nachalo pravogo plecha
                    this.add(new helix(ls1, ls2, rs1, rs2, seq));
                    ls1 = ls2 = rs1 = rs2 = -1;
                    counter = 0;
                  }
                  if (match[i - k + t][j + k - t] == 1)
                  {
                    counter--;
                    continue;
                  }
                }
                else break;
              }
            }
          }
          counter = 0;
        }
      }
    }
//=============================================
    // udalyaem shpilki sami k sebe, t.e. vida (3-5:5-3 ili 3-5:6-4)
    for (int i = 0; i < this.size(); i++)
    {
      helix sh = (helix) this.get(i);
      if ((sh.A <= sh.C) && (sh.C <= sh.B))
      {
        this.remove(i--);
      }
    }

    // udalyaem odinakovye shpilki
    for (int i = 0; i < this.size(); i++)
    {
      helix sh = (helix) this.get(i);
      for (int j = i + 1; j < this.size(); j++)
      {
        helix sh1 = (helix) this.get(j);
        if ( (sh.A == sh1.C) && (sh.B == sh1.D) && (sh.C == sh1.A) &&
            (sh.D == sh1.B))
       {
          this.remove(j--);
        }
      }
    }

    // shpilki bez petli ili s petlei < 3bp, vida 3-5:8-6
    // nado smotret' mozhno li rasplesti nemnogo spiral', chtoby sdelat' normal'nuyu shpil'ku
    // esli nel'zya, to udalyaem
    for (int i = 0; i < this.size(); i++)
    {
      helix sh = (helix) this.get(i);
      if (sh.C - sh.B <= minloop)
      {
        int helix_length = sh.Left.length();
        if (helix_length == minhelix) // esli dlina spirali = minimal'noi dline, to nel'zya rasplesti
          this.remove(i--);
        else // esli dlina spirali bol'she minimal'noi i mozhno poprobovat' rasplesti
          for (int j=1; j<helix_length; j++)
          {
            helix_length--;
            if (helix_length >= minhelix)  // t.e. mozhno rasplesti
            {
              if (match[sh.B - j][sh.C + j] == 2)  // esli teper' spiral' zakanchivaetsya na kanonicheskuyu paru
              {
                sh.B = sh.B-j;
                sh.C = sh.C+j;
                sh.Left = seq.substring(sh.A, sh.B + 1); // levoe plecho
                sh.Right = sh.reverse(seq.substring(sh.C, sh.D + 1)); // perevernutoe pravoe
                sh.energy = sh.energy_calculation(sh.Left, sh.Right);
                break;
              }
              else continue;  // esli zakanchivaetsya na (g,u)
            }
            else  this.remove(i--);
          }
      }
    }
    // udalyaem malozhivushe shpilki
    //added by A.Favorov
 /*   for (int i = 0; i < this.size(); i++)
    {
      helix sh = (helix)this.get(i);
      if(sh.k_break > (1./min_helix_living_time))
        this.remove(i--);
    }//*/

    // udalyaem spirali s energiei bolshe min_energy
     for (int i = 0; i < this.size(); i++)
    {
      helix sh = (helix)this.get(i);
      if(sh.energy > (min_energy))
        this.remove(i--);
    }


    //System.out.println("");
    sort(); // uporyadochivaem vector helices
    heliset_size = this.size();

    //=====================================
    // ischem spirali, kotorye mogut obrazovat' structuru so skol'zhyaschei petlei
    // i kladem ih v slid_heliset
    slid_heliset_process = new Vector();
    slid_heliset_process.add(0, null); // kladem pustoi element, chtoby indeksy u sk. spiralei byli >0
    // helix1 vyshe helix2
    for (int i = 0; i < heliset_size; i++)
    {
      helix helix1 = (helix) get(i);
      int A1 = helix1.A; // nachalo levogo plecha
      int B1 = helix1.B; // konec levogo plecha
      int C1 = helix1.C; // nachalo pravogo plecha
      int D1 = helix1.D; // konec pravogo plecha
      for (int j = 0; j < heliset_size; j++)
      {
        helix helix2 = (helix) get(j);
        int A2 = helix2.A; // nachalo levogo plecha
        int B2 = helix2.B; // konec levogo plecha
        int C2 = helix2.C; // nachalo pravogo plecha
        int D2 = helix2.D; // konec pravogo plecha
        // spirali idut odna nad drugoi (peresekayutsya libo levye plechi, libo pravye)
        if ((A2 - A1 < D1 - D2) && // uslovia naklona vlevo
            (D1 >= C2) && (C2 > C1) && // pravye plechi peresekayutsya
 //           (B2 < A1) && // a levye ne peresekayutsya
            (C2 - C1 >= minpartialhelix) && // bez peresecheniya ostaetsya sparennymi ne men'she (minhelix-1) par
            (D2 - D1 >= minpartialhelix))
        {
          slid_heliset_process.add(new slid_helix(i, j, helix1, helix2, false));
          continue;
        }
        if (( A2 - A1 > D1 - D2) && // naklon vpravo
            (B1 > B2) && (B2 >= A1) &&// levye plechi peresekayutsya
 //           (C2 > D1) && // a pravye net
            (B1 - B2 >= minpartialhelix) &&
            (A1 - A2 >= minpartialhelix))
        {
          slid_heliset_process.add(new slid_helix(i, j, helix1, helix2, false));
          continue;
        }

        // spirali idut posledovatel'no (peresekayutsya raznoimennye plechi)
        if (( A2 > C1) && (A2 <= D1) && (B2 > D1) &&// levoe plecho sp. i peresekaetsya s pravym sp. j
            (A2 - C1 >= minpartialhelix) &&
            (B2 - D1 >= minpartialhelix))
        {
          slid_heliset_process.add(new slid_helix(i, j, helix1, helix2, true));
          continue;
        }
      }
    }

    //================================
    // ubiraem odinakovye iz slid_heliset
  /*  for (int i = 1; i < slid_heliset.size(); i++)
    {
      slid_helix sh = (slid_helix) slid_heliset.get(i);
      for (int j = i + 1; j < slid_heliset.size(); j++)
      {
        slid_helix sh1 = (slid_helix) slid_heliset.get(j);
        if ( sh.slip_parts[0] == sh1.slip_parts[1] && sh.slip_parts[1] == sh1.slip_parts[0])
        {
          slid_heliset.remove(j--);
        }
      }
    }//*/
    //=================================
    // zapolnyaem massiv Compatible
    Compatible = new int[heliset_size][heliset_size];
    // dlya pravil'nyh spiralei
    for (int i = 0; i < this.heliset_size; i++)
    {
      helix sh = (helix) this.get(i);
      for (int j = i + 1; j < heliset_size; j++)
      {
        helix sh1 = (helix) this.get(j);
        Compatible[i][j] = Compatible[j][i] = isCompatible(sh, sh1);
      }
    }
    // dlya skol'zyaschih petel'
    // delaem iz slid_heliset_process massiv slid_heliset
    // i zapolnyaem Compatible
    slid_heliset = new slid_helix[slid_heliset_process.size()];
    slid_heliset = (slid_helix[]) slid_heliset_process.toArray(slid_heliset);
    for (int i = 1; i < slid_heliset_process.size(); i++)
    {
      slid_helix sh = (slid_helix) slid_heliset_process.get(i);
      //slid_heliset[i] = sh;
      Compatible[sh.slip_parts[0]][sh.slip_parts[1]] =
          Compatible[sh.slip_parts[1]][sh.slip_parts[0]] = i;
    }

    slid_heliset_process = null;

    //======================================
    // ischem svyazannye pary spiralei
    bind_pairs = new int[heliset_size][heliset_size];
    for (int i = 0; i < heliset_size; i++)
      for (int j = 0; j < heliset_size; j++)
      {
        bind_pairs[i][j] = -1;
      }
    for (int i = 0; i < heliset_size; i++)
      for (int j = i + 1; j < heliset_size; j++)
      { // esli spirali i i j sovmestimy
        if (Compatible[i][j] == 0)
        {
          // to ischem spiral', s kotoroi oni obe nesovmestimy
          for (int k = 0; k < heliset_size; k++)
            if (k != i && k != j &&
                Compatible[i][k] == -1 &&
                Compatible[j][k] == -1)
            {
              helix sh_i = (helix)this.get(i);
              helix sh_j = (helix)this.get(j);
              helix sh_k = (helix)this.get(k);
              // proveryaem kak nesovmestimy
              // plechi peresekayutsya sverhu
              if ( (sh_k.A <= sh_i.A && sh_i.A <= sh_k.B && // peresechenie levyh plechei k i i
                    sh_i.D < sh_k.C && // pravye ne peresekayutsya
                    sh_k.D >= sh_j.D && sh_j.D >= sh_k.C && // peresechenie pravyh plechei k i j
                    sh_j.A > sh_k.B) || // levye ne peresekayutsya
                  // ili naoborot
                  (sh_k.A <= sh_j.A && sh_j.A <= sh_k.B && // peresechenie levyh plechei k i j
                   sh_j.D < sh_k.C && // pravye ne peresekayutsya
                   sh_k.D >= sh_i.D && sh_i.D >= sh_k.C && // peresechenie pravyh plechei k i i
                   sh_i.A > sh_k.B))
                bind_pairs[i][j] = bind_pairs[j][i] = k;
                // plechi peresekayutsya snuzu
              if ( (sh_k.A <= sh_i.D && sh_i.D <= sh_k.B && // peresechenie levogo plecha k i pravogo i
                    //sh_i.A < sh_k.C && // ?
                    sh_k.D >= sh_j.A && sh_k.D <= sh_j.B) || // peresechenie pravogo plecha k i levogo j
                  //  sh_j.A > sh_k.B) || // ?
                  // ili naoborot
                  (sh_k.A <= sh_j.D && sh_j.D <= sh_k.B && // peresechenie levogo plecha k i pravogo j
                   sh_k.D >= sh_i.A && sh_k.D <= sh_i.B)) // peresechenie pravogo plechf k i levogo i

                bind_pairs[i][j] = bind_pairs[j][i] = k; //*/
            }

        }
     }


  }

//===============================================
  // proverka shpilek na sovmestimost
  public int isCompatible(helix a, helix b)
  {
    if ((a.B <= b.A) && (b.A <= a.D) && (a.C <= b.D)) return -1;
    if ((b.B <= a.A) && (a.A <= b.C) && (b.C <= a.D)) return -1;
    if (((a.A <= b.A) && (b.A <= a.B)) || ((b.A <= a.A) && (a.A <= b.B)))
       return -1;
    if (((b.C <= a.D) && (a.D <= b.D)) || ((a.C <= b.D) && (b.D <= a.D)))
      return -1;
    return 0;
  }


  //============================================
  // govorit est' li danna spiral' v heliset i vozvraschaet nomer
  int contain(helix a)
  {
    int A1 = a.A;
    int B1 = a.B;
    int C1 = a.C;
    int D1 = a.D;
    for (int i = 0; i< this.size(); i++)
    {// perebiraem heliset
      // vytaskivaem spiral' iz heliset
      helix b = (helix) get(i);
      int A2 = b.A;
      int B2 = b.B;
      int C2 = b.C;
      int D2 = b.D;
      if (A1 == A2 && B1 == B2 && C1 == C2 && D1 == D2) // nashli takuyu zhe spiral'
        return i;
      // esli a yavlyaetsya podspiral'yu b
      if (A1 == A2 && B1 <= B2 && C1 >= C2 && D1 == D2)
        return i;
      if (A1 >= A2 && B1 == B2 && C1 == C2 && D1 <= D2)
        return i;
      if (A1 >= A2 && D1 <= D2 && A1-A2 == D2-D1 &&
          B1 <= B2 && C1 >= C2 && B2-B1 == C1-C2)
        return i;
    }
    return -1;
  }

//====================================
  // funkciya, kotoraya po nomeru spirali smotrit s kem ona obrazuet sk. petlyu,
  // s kakih koncov i plechei i vydaet massiv, sostoyaschii iz dlina spirali bez sk. petel'
  // i energiyu etogo kusochka
  // 0-oi element: dlina
  // 1-yi: energiya
  double[] get_length_energi_without_slid(int t, boolean[] exists)
  {
    double result[] = new double[2]; // 0-oi element: dlina; 1-yi: energiya
    helix helix_t = (helix) get(t); // hochet obrazovat'sya

  // pomnyat nomera slid spiralei, obrazuyuschoh s t-oi sk. petli s 4-h raznyh koncov
    int top_right = 0, top_left = 0, bottom_right = 0, bottom_left = 0;

    double split_energy = helix_t.energy;
    int n = helix_t.Left.length();
    String t_left = helix_t.Left;
    String t_right = helix_t.Right;

    //===================================================
    //probegaem vsie spirali, ishem sobratiev po skolzishim koncam i zapominaem ih
    for (int j = 0; j < heliset_size; j++)
    {
      int Compat = Compatible[t][j];
      if ((exists[j]) && (Compat > 0))
      {
        slid_helix slid_helix1 = (slid_helix) slid_heliset[Compat];
        helix helix_j = (helix) get(j);
        //========================================================
        // esli peresekayutsya raznoimennye plechi
        if (slid_helix1.bottom_slid)
        {
          // pravoe plecho t-oi i levoe j-oi
          if (helix_t.A < helix_j.A)
          {
            bottom_left = Compat;
            continue;
          }
          else// levoe plecho t-oi i pravoe j-oi
          {
              bottom_right = Compat;
               continue;
          }
        }
        //==========================================================
        else // esli peresekayutsya odnoimennye plechi
        // esli sk. petlya na pravom pleche (left_slid = false), t.e. peresechenie levyh
        if (!slid_helix1.left_slid)
        {
        // i j-aya spiral' vyshe t-oi
          if (helix_t.A < helix_j.A)
          {
            top_left = Compat;
            continue;
          }
          // esli peresechenie levyh i j-aya spiral' nizhe t-oi
          else // helix_t.A > helix_j.A)
          {
            bottom_left = Compat;
            continue;
          }
        }
        // esli sk. petlya na levom pleche (left_slid = true), t.e. peresechenie pravyh plechei
        else
        {
          // i j-aya spiral' vyshe t-oi
          if (helix_t.A < helix_j.A)
          {
            top_right = Compat;
            continue;
          }
          // esli peresechenie pravyh i j-aya spiral' nizhe t-oi
          else // helix_t.A > helix_j.A)
          {
            bottom_right = Compat;
            continue;
          }
        }
      }
    }
    //===========================================
    // esli 2 sk. petli, to vychitaem bol'shuyu dlinu i bol'shuyu energiyu
    // esli sverhu
    if (top_right > 0 && top_left > 0)
    {
      // pravaya petlya
      slid_helix right = (slid_helix) slid_heliset[top_right];
      // levaya petlya
      slid_helix left = (slid_helix) slid_heliset[top_left];
      // esli na pravom pleche peresechenie bol'she, to vychitaem peresechenie iz dliny i energii
      if (right.cross > left.cross)
      {
        split_energy -= right.energy_cross;
        n -= right.cross;
        t_left = t_left.substring(0, t_left.length() - right.cross + 1);
        t_right = t_right.substring(0, t_right.length() - right.cross + 1);
      }
      else // esli na levom bol'she, to vychitaem peresechenie iz dliny i energii
      {
        split_energy -= left.energy_cross;
        n -= left.cross;
        t_left = t_left.substring(0, t_left.length() - left.cross + 1);
        t_right = t_right.substring(0, t_right.length() - left.cross + 1);
      }
    }
    else // kogda petlya odna ili niodnoi
    {
      if (top_right > 0)
      {
        slid_helix right = (slid_helix) slid_heliset[top_right];
        split_energy -= right.energy_cross;
        n -= right.cross;
        t_left = t_left.substring(0, t_left.length() - right.cross + 1);
        t_right = t_right.substring(0, t_right.length() - right.cross + 1);
      }
      if (top_left > 0)
      {
        slid_helix left = (slid_helix) slid_heliset[top_left];
        split_energy -= left.energy_cross;
        n -= left.cross;
        t_left = t_left.substring(0, t_left.length() - left.cross + 1);
        t_right = t_right.substring(0, t_right.length() - left.cross + 1);
      }
    }
    // esli snizu
    if (bottom_right > 0 && bottom_left > 0)
    {
      // pravaya petlya
      slid_helix right = (slid_helix) slid_heliset[bottom_right];
      // levaya petlya
      slid_helix left = (slid_helix) slid_heliset[bottom_left];
      // esli na pravom pleche peresechenie bol'she, to vychitaem peresechenie iz dliny i energii
      if (right.cross > left.cross)
      {
        split_energy -= right.energy_cross;
        n -= right.cross;
        t_left = t_left.substring(right.cross);
        t_right = t_right.substring(right.cross);
      }
      else // esli na levom bol'she,to vychitaem peresechenie iz dliny i energii
      {
        split_energy -= left.energy_cross;
        n -= left.cross;
        t_left = t_left.substring(left.cross);
        t_right = t_right.substring(left.cross);
      }
    }
    else // kogda petlya odna ili niodnoi
    {
      if (bottom_right > 0)
      {
        slid_helix right = (slid_helix) slid_heliset[bottom_right];
        split_energy -= right.energy_cross;
        n -= right.cross;
        t_left = t_left.substring(right.cross);
        t_right = t_right.substring(right.cross);
      }
      if (bottom_left > 0)
      {
        slid_helix left = (slid_helix) slid_heliset[bottom_left];
        split_energy -= left.energy_cross;
        n -= left.cross;
        t_left = t_left.substring(left.cross);
        t_right = t_right.substring(left.cross);
      }
    }
    split_energy = helix_t.energy_calculation(t_left,t_right);
    result[0] = n; // zapominaem dlinu
    result[1] = split_energy; // i energiyu
    return result;
  }


//====================================
  // funkciya, kotoraya po nomeru spirali smotrit s kem ona obrazuet sk. petlyu,
  // s kakih koncov i plechei, a takzhe smotrit ne obrazuyut li
    //  te spirali v svoyu ochered' eshe s kem-to, i esli posle vycheta peresecheniy
    // ih dlina budet <= 0, to takyuy spiral' nel'zya obrazovyvat'
    // vozvrzschaet dlinu spirali bez sk. petli
    // vozvraschaet 0, esli etu spiral' nel'zya obrazovat'

    int get_length_without_slid(int t, boolean[] exists)
     {
       exists[t] = true; // pokazyvaem, chto t-ya hochet obrazovat'sya
       // i smotrim kak bubut sebya vesti ostal'nye spirali
       helix helix_t = (helix) get(t); // hochet obrazovat'sya
       int n = helix_t.Left.length(); // dlina tekuschei spirali
       // pomnyat nomera slid spiralei (iz slid_heliset), obrazuyuschoh s t-oi sk. petli s 4-h raznyh koncov
       int top_right = 0, top_left = 0, bottom_right = 0, bottom_left = 0;
       // pomnyat nomera spiralei (iz heliset)
       int top_right_helix = -1, top_left_helix = -1,
           bottom_right_helix = -1, bottom_left_helix = -1;
       // dliny spiralei, s kotorymi t-aya obrazuet sk. petli
       int n_top_right = 0, n_top_left = 0, n_bottom_right = 0, n_bottom_left = 0;
       //===================================================
       //probegaem vsie spirali, ishem sobratiev po skolzishim koncam i zapominaem ih
       for (int j = 0; j < heliset_size; j++)
       {
         int Compat = Compatible[t][j];
         if ((exists[j]) && (Compat > 0))
         {
           slid_helix slid_helix1 = (slid_helix) slid_heliset[Compat];
           helix helix_j = (helix) get(j);
           //=================================================
           // esli peresekayutsya raznoimennye plechi
           if (slid_helix1.bottom_slid)
           {
             // pravoe plecho t-oi i levoe j-oi
             if (helix_t.A < helix_j.A)
             {
               bottom_left = Compat;
               bottom_left_helix = j;
               continue;
             }
             else // levoe plecho t-oi i pravoe j-oi
             {
               bottom_right = Compat;
               bottom_right_helix = j;
               continue;
             }
           }
           //============================================
           // esli sk. petlya na pravom pleche (left_slid = false), t.e. peresechenie levyh plechei
           if (!slid_helix1.left_slid)
           {
             // i j-aya spiral' vyshe t-oi
             if (helix_t.A < helix_j.A)
             {
               top_left = Compat;
               top_left_helix = j;
               continue;
             }
             // esli peresechenie levyh i j-aya spiral' nizhe t-oi
             else // helix_t.A > helix_j.A)
             {
               bottom_left = Compat;
               bottom_left_helix = j;
               continue;
             }
           }
           // esli sk. petlya na levom pleche (left_slid = true), t.e. peresechenie pravyh plechei
           else
           {
             // i j-aya spiral' vyshe t-oi
             if (helix_t.A < helix_j.A)
             {
               top_right = Compat;
               top_right_helix = j;
               continue;
             }
             // esli peresechenie pravyh i j-aya spiral' nizhe t-oi
             else // helix_t.A > helix_j.A)
             {
               bottom_right = Compat;
               bottom_right_helix = j;
               continue;
             }
           }
         }
      }
      exists[t] = false; // delaem kak bylo

       //===========================================
       // esli 2 sk. petli, to vychitaem bol'shuyu dlinu i bol'shuyu energiyu
       // esli sverhu
       // massiv dlin, kotorye nado vychest' iz top_right_helix, top_left_helix,
       // bottom_right_helix i bottom_left_helix s 4 koncov
       // 0 - bottom_left, 1 - bottom_right, 2 - top_left, 3 - top_right
       int[] top_right_lengths = new int[4];
       int[] top_left_lengths = new int[4];
       int[] bottom_right_lengths = new int[4];
       int[] bottom_left_lengths = new int[4];
       for (int j = 0; j < 4; j++)
       {
         top_right_lengths[j] = 0;
         top_left_lengths[j] = 0;
         bottom_right_lengths[j] = 0;
         bottom_left_lengths[j] = 0;
       }
       if (top_right > 0 && top_left > 0)
       {
         // pravaya petlya
         slid_helix right = (slid_helix) slid_heliset[top_right];
         // levaya petlya
         slid_helix left = (slid_helix) slid_heliset[top_left];
         // dlina pravoi spirali
         n_top_right = ((helix) get(top_right_helix)).Left.length();
         // spiral' levaya
         n_top_left = ((helix) get(top_left_helix)).Left.length();

         top_right_lengths[1] = right.cross;
         top_left_lengths[0] = left.cross;
         // ischem sk. petli dlya pravoi i levoi spiralei
         // esli nahodim, vychitaem
         for (int j = 0; j < heliset_size; j++) //probegaem vsie, ishem sobratiev po skolzishim koncam
         {
           int Compat_right = Compatible[top_right_helix][j];
           int Compat_left = Compatible[top_left_helix][j];
           if ( (exists[j]) && (Compat_right > 0))
           {
             // nashli sk. petlyu dlya spirali top_right
             slid_helix right1 = (slid_helix) slid_heliset[Compat_right];
             // smortim peresekayutsya li oni neodnoimennymi plechami
             // esli da, sk. petli poluchayutsya na raznyh plechah snizhu u top_right
             // i vychitat' nado bolshee iz peresechenii
             if (right1.bottom_slid)
             {
               // peresechenie vnizu sleva
               top_right_lengths[0] = right1.cross;
 /*              if (right1.cross > right.cross)
                 n_top_right -= right1.cross;// otnimaem peresechenie s right1
               else n_top_right -= right.cross; // inache otnimaem dlinu peresecheniya s t-oi*/
             }
             else // znachit petlya sverhu u top_right
             {
               if (right1.left_slid) // peresechenie sprava
                 top_right_lengths[3] = right1.cross;
               else // peresechenie sleva
                 top_right_lengths[2] = right1.cross;
//             n_top_right -= right1.cross; // otnimaem peresechenie s right1
//             n_top_right -= right.cross; // i otnimaem dlinu peresecheniya s t-oi
             }
             continue;
           }
           if ( (exists[j]) && (Compat_left > 0))
           {
             // nashli sk. petlyu dlya spirali top_left
             slid_helix left1 = (slid_helix) (slid_heliset[Compat_left]);
             // smortim peresekayutsya li oni neodnoimennymi plechami
             // esli da, sk. petli poluchayutsya na raznyh plechah snizhu u top_left
             // i vychitat' nado bolshee iz peresechenii
             if (left1.bottom_slid)
             {
               // peresechenie vnizu sprava
               top_left_lengths[0] = left1.cross;
/*               if (left1.cross > left.cross)
                 n_top_left -= left1.cross;// otnimaem peresechcnie s left1
               else n_top_left -= left.cross;// inache otnimaem dlinu peresecheniya s t-oi*/
             }
             else // znachit petlya sverhu u top_left
             {
               if (left1.left_slid) // peresechenie sprava
                 top_left_lengths[3] = left1.cross;
               else
                 top_left_lengths[2] = left1.cross;
//               n_top_left -= left1.cross; //otnimaem peresechcnie s left1
//               n_top_left -= left.cross; // i otnimaem dlinu peresecheniya s t-oi
             }
             continue;
           }
         }
         // korrektiruem dlinu t-oi spirali sverhu
         // esli na pravom pleche peresechenie bol'she, to vychitaem peresechenie iz dliny
         if (right.cross > left.cross)
         {
           n -= right.cross;
         }
         else // esli na levom bol'she, to vychitaem peresechenie iz dliny
         {
           n -= left.cross;
         }
       }
       else // kogda petlya odna ili niodnoi
       {
         if (top_right > 0)
         {
           // dlina pravoi spirali
           n_top_right = ( (helix) get(top_right_helix)).Left.length();
           // peresechenie s t-oi
           int cross = ( (slid_helix) (slid_heliset[top_right])).cross;
           top_right_lengths[1] = cross;
         // ischem sk. petli dlya pravoi i levoi spiralei
           // esli nahodim, vychitaem
           for (int j = 0; j < heliset_size; j++) //probegaem vsie, ishem sobratiev po skolzishim koncam
           {
             int Compat_right = Compatible[top_right_helix][j];
             if ( (exists[j]) && (Compat_right > 0))
             {
             // nashli sk. petlyu dlya spirali top_right
             slid_helix right1 = (slid_helix) slid_heliset[Compat_right];
             // smortim peresekayutsya li oni neodnoimennymi plechami
             // esli da, sk. petli poluchayutsya na raznyh plechah snizhu u top_right
             // i vychitat' nado bolshee iz peresechenii
             if (right1.bottom_slid)
             {
               // peresechenie vnizu sleva
               top_right_lengths[0] = right1.cross;
 /*              if (right1.cross > right.cross)
                 n_top_right -= right1.cross;// otnimaem peresechenie s right1
               else n_top_right -= right.cross; // inache otnimaem dlinu peresecheniya s t-oi*/
             }
             else // znachit petlya sverhu u top_right
             {
               if (right1.left_slid) // peresechenie sprava
                 top_right_lengths[3] = right1.cross;
               else // peresechenie sleva
                 top_right_lengths[2] = right1.cross;
//             n_top_right -= right1.cross; // otnimaem peresechenie s right1
//             n_top_right -= right.cross; // i otnimaem dlinu peresecheniya s t-oi
             }
             continue;
             }
           }
           // korrektiruem dlinu t-oi
           n -= cross;
         }
         if (top_left > 0)
         {
           // spiral' levaya
           n_top_left = ( (helix) get(top_left_helix)).Left.length();
           // peresechenie s t-oi
           int cross = ( (slid_helix) (slid_heliset[top_left])).cross;
           top_left_lengths[0] = cross;
           // ischem sk. petli dlya pravoi i levoi spiralei
           // esli nahodim, vychitaem
           for (int j = 0; j < heliset_size; j++) //probegaem vsie, ishem sobratiev po skolzishim koncam
           {
             int Compat_left = Compatible[top_left_helix][j];
             if ( (exists[j]) && (Compat_left > 0))
             {
               // nashli sk. petlyu dlya spirali top_left
               slid_helix left1 = (slid_helix) (slid_heliset[Compat_left]);
               // smortim peresekayutsya li oni neodnoimennymi plechami
               // esli da, sk. petli poluchayutsya na raznyh plechah snizhu u top_left
               // i vychitat' nado bolshee iz peresechenii
               if (left1.bottom_slid)
               {
                 // peresechenie vnizu sprava
                 top_left_lengths[0] = left1.cross;
/*               if (left1.cross > left.cross)
                  n_top_left -= left1.cross;// otnimaem peresechcnie s left1
                  else n_top_left -= left.cross;// inache otnimaem dlinu peresecheniya s t-oi*/
               }
               else // znachit petlya sverhu u top_left
               {
                 if (left1.left_slid) // peresechenie sprava
                   top_left_lengths[3] = left1.cross;
                 else
                   top_left_lengths[2] = left1.cross;
//               n_top_left -= left1.cross; //otnimaem peresechcnie s left1
//               n_top_left -= left.cross; // i otnimaem dlinu peresecheniya s t-oi
               }
               continue;
             }
           }
           // korrektiruem dlinu t-oi
           n -= cross;
         }
       }
       // esli snizu
       if (bottom_right > 0 && bottom_left > 0)
       {
         // pravaya petlya
         slid_helix right = (slid_helix) slid_heliset[bottom_right];
         // levaya petlya
         slid_helix left = (slid_helix) slid_heliset[bottom_left];
         // dlina pravoi spirali
         n_bottom_right = ((helix) get(bottom_right_helix)).Left.length();
         // dlina levoi spirali
         n_bottom_left = ((helix) get(bottom_left_helix)).Left.length();
         if (right.bottom_slid) bottom_right_lengths[0] = right.cross;
         else bottom_right_lengths[3] = right.cross;
         if (left.bottom_slid) bottom_left_lengths[1] = left.cross;
         else bottom_left_lengths[2] = left.cross;
         // ischem sk. petli dlya pravoi i levoi spiralei
         // esli nahodim, vychitaem
         for (int j = 0; j < heliset_size; j++) //probegaem vsie, ishem sobratiev po skolzishim koncam
         {
           int Compat_right = Compatible[bottom_right_helix][j];
           int Compat_left = Compatible[bottom_left_helix][j];
           if ( (exists[j]) && (Compat_right > 0))
           {
             // nashli sk. petlyu dlya spirali bottom_right
             slid_helix right1 = (slid_helix) slid_heliset[Compat_right];
             // smortim peresekayutsya li oni neodnoimennymi plechami
             // esli da, sk. petli poluchayutsya na raznyh plechah snizhu u bottom_right
             // i vychitat' nado bolshee iz peresechenii
             if (right1.bottom_slid)
             {
               bottom_right_lengths[1] = right1.cross;
             }
             else // zavisit ot otnosheniya s t-oi
             {
               if (right.bottom_slid && right1.left_slid) // peresechenie sverhu sprava
                 bottom_right_lengths[3] = right1.cross;
               if (right.bottom_slid && !right1.left_slid) // peresechenie sverhu sleva
                 bottom_right_lengths[2] = right1.cross;
               if (!right.bottom_slid && right1.left_slid) // snizu sprava
                 bottom_right_lengths[1] = right1.cross;
               if (!right.bottom_slid && !right1.left_slid) // snizhu sleva
                 bottom_right_lengths[0] = right1.cross;

//               n_bottom_right -= right1.cross; // otnimaem peresecheni s right1
//               n_bottom_right -= right.cross; // i peresechenie s t-oi
             }
             continue;
           }
           if ( (exists[j]) && (Compat_left > 0))
           {
             slid_helix left1 = (slid_helix) slid_heliset[Compat_left];
             if (left1.bottom_slid)
             {
               bottom_left_lengths[0] = left1.cross;
             }
             else // zavisit ot otnosheniya s t-oi
             {
               if (left.bottom_slid && left1.left_slid) // sverhu sprava
                 bottom_left_lengths[3] = left1.cross;
               if (left.bottom_slid && !left1.left_slid) // sverhu sleva
                 bottom_left_lengths[2] = left1.cross;
               if (!left.bottom_slid && left1.left_slid) // snizu sprava
                 bottom_left_lengths[1] = left1.cross;
               if (!left.bottom_slid && !left1.left_slid) // snizu sleva
                 bottom_left_lengths[0] = left1.cross;
//               n_bottom_left -= left1.cross; //otnimaem peresechenie s left1
//               n_bottom_left -= left.cross; // i s t-oi
             }
             continue;
           }
         }
         // korrektiruem dline t-oi snizu
       // esli na pravom pleche peresechenie bol'she, to vychitaem peresechenie iz dliny i energii
         if (right.cross > left.cross)
         {
           n -= right.cross;
         }
         else // esli na levom bol'she,to vychitaem peresechenie iz dliny i energii
         {
           n -= left.cross;
         }
       }
       else // kogda petlya odna ili niodnoi
       {
         if (bottom_right > 0)
         {
           // dlina pravoi spirali
           n_bottom_right = ( (helix) get(bottom_right_helix)).Left.length();
           slid_helix right = (slid_helix) slid_heliset[bottom_right];
           int cross = right.cross;
           if (right.bottom_slid) bottom_right_lengths[0] = cross;
           else bottom_right_lengths[3] = cross;

           for (int j = 0; j < heliset_size; j++) //probegaem vsie, ishem sobratiev po skolzishim koncam
           {
             int Compat_right = Compatible[bottom_right_helix][j];
             if ( (exists[j]) && (Compat_right > 0))
             {
               // nashli sk. petlyu dlya spirali bottom_right
               slid_helix right1 = (slid_helix) slid_heliset[Compat_right];
               // smortim peresekayutsya li oni neodnoimennymi plechami
               // esli da, sk. petli poluchayutsya na raznyh plechah snizhu u bottom_right
               // i vychitat' nado bolshee iz peresechenii
               if (right1.bottom_slid)
               {
                 bottom_right_lengths[1] = right1.cross;
               }
               else // zavisit ot otnosheniya s t-oi
               {
                 if (right.bottom_slid && right1.left_slid) // peresechenie sverhu sprava
                   bottom_right_lengths[3] = right1.cross;
                 if (right.bottom_slid && !right1.left_slid) // peresechenie sverhu sleva
                   bottom_right_lengths[2] = right1.cross;
                 if (!right.bottom_slid && right1.left_slid) // snizu sprava
                   bottom_right_lengths[1] = right1.cross;
                 if (!right.bottom_slid && !right1.left_slid) // snizhu sleva
                   bottom_right_lengths[0] = right1.cross;

   //               n_bottom_right -= right1.cross; // otnimaem peresecheni s right1
   //               n_bottom_right -= right.cross; // i peresechenie s t-oi
               }
               continue;
             }
           }
           n -= cross;
         }
         if (bottom_left > 0)
         {
           // dlina levoi spirali
           n_bottom_left = ( (helix) get(bottom_left_helix)).Left.length();
           slid_helix left = (slid_helix) slid_heliset[bottom_left];
           int cross = left.cross;
           if (left.bottom_slid) bottom_left_lengths[1] = left.cross;
           else bottom_left_lengths[2] = left.cross;

           for (int j = 0; j < heliset_size; j++) //probegaem vsie, ishem sobratiev po skolzishim koncam
           {
             int Compat_left = Compatible[bottom_left_helix][j];
             if ( (exists[j]) && (Compat_left > 0))
             {
               slid_helix left1 = (slid_helix) slid_heliset[Compat_left];
             if (left1.bottom_slid)
             {
               bottom_left_lengths[0] = left1.cross;
             }
             else // zavisit ot otnosheniya s t-oi
             {
               if (left.bottom_slid && left1.left_slid) // sverhu sprava
                 bottom_left_lengths[3] = left1.cross;
               if (left.bottom_slid && !left1.left_slid) // sverhu sleva
                 bottom_left_lengths[2] = left1.cross;
               if (!left.bottom_slid && left1.left_slid) // snizu sprava
                 bottom_left_lengths[1] = left1.cross;
               if (!left.bottom_slid && !left1.left_slid) // snizu sleva
                 bottom_left_lengths[0] = left1.cross;
//               n_bottom_left -= left1.cross; //otnimaem peresechenie s left1
//               n_bottom_left -= left.cross; // i s t-oi
             }

               continue;
             }
           }
           n -= cross;
         }
       }
       // korrektiruem dlinu n_top_right
       if (top_right_lengths[0] > top_right_lengths[1])
         n_top_right -= top_right_lengths[0];
       else
         n_top_right -= top_right_lengths[1];
       if (top_right_lengths[2] > top_right_lengths[3])
         n_top_right -= top_right_lengths[2];
       else
         n_top_right -= top_right_lengths[3];
         // korrektiruem dlinu n_top_left
       if (top_left_lengths[0] > top_left_lengths[1])
         n_top_left -= top_left_lengths[0];
       else
         n_top_left -= top_left_lengths[1];
       if (top_left_lengths[2] > top_left_lengths[3])
         n_top_left -= top_left_lengths[2];
       else
         n_top_left -= top_left_lengths[3];
         // korrektiruem dlinu n_bottom_right
       if (bottom_right_lengths[0] > bottom_right_lengths[1])
         n_bottom_right -= bottom_right_lengths[0];
       else n_bottom_right -= bottom_right_lengths[1];
       if (bottom_right_lengths[2] > bottom_right_lengths[3])
        n_bottom_right -= bottom_right_lengths[2];
      else n_bottom_right -= bottom_right_lengths[3];
      // korrektiruem dlinu n_bottom_left
       if (bottom_left_lengths[0] > bottom_left_lengths[1])
         n_bottom_left -= bottom_left_lengths[0];
       else n_bottom_left -= bottom_left_lengths[1];
       if (bottom_left_lengths[2] > bottom_left_lengths[3])
        n_bottom_left-= bottom_left_lengths[2];
      else n_bottom_left -= bottom_left_lengths[3];

       // proveryaem, chto esli obrazuetsya eta spiral', to nikakaya s nei smezhnaya ne budet dlina <= 0
       // esli budet, to vozvraschaem 0
       if (n_top_right < 0 || n_top_left < 0 ||
           n_bottom_right < 0 || n_bottom_left < 0) n = 0;
       return n;
     }


//===================================
  // uporyadochivaem vector helices po koncu spirali rs2
  void sort()
  {
    for (int i = 0; i < this.size(); i++)
    {
      helix sh = (helix) this.get(i);
      for (int j = i + 1; j < this.size(); j++)
      {
        helix sh1 = (helix)this.get(j);
        if (sh.D > sh1.D)
        {
          this.remove(j);
          this.add(i, sh1);
        }
      }
    }
    int t = -1;
    for (int i = 0; i < this.size(); i++)
    {
      helix sh = (helix)this.get(i);
      int ls = sh.D;
      for (int j = i + 1; j < this.size(); j++)
      {
        helix sh1 = (helix)this.get(j);
        if (sh1.D < ls)
        {
          ls = sh1.D;
          t = j;
        }
      }
      if (t != -1)
      {
        helix sh1 = (helix)this.get(t);
        this.remove(t);
        this.add(i, sh1);
        t = -1;
      }
      else continue;
    }
  }

   //===================================================
   // zapolnyaem grafik dlya spiralei
   void add_to_graph (double CurrTime, double LiveTime, Vector exist_helises)
   {
     // ischem mesto v massive graph
     double time_start = CurrTime; //nachalo zhizni struktury
     double time_end = CurrTime + LiveTime;
     if (CurrTime + LiveTime > rna1.Teta) time_end = rna1.Teta; // konec zhizni
     int n_start = (int)((rna1.graph_step*time_start)/rna1.Teta); // element massiva graph, kuda popadaet nachalo
     int n_end = (int)((rna1.graph_step*time_end)/rna1.Teta); // element massiva graph, kuda popadaet konec
     for (int i = n_start; i<= n_end; i++)
     {
       double next_time = 0;
       if (i+1 < rna1.graph_step)
       {
         next_time = ( (i + 1) * rna1.Teta) / rna1.graph_step; // vremya, sootvetstvuyuschee elementu i+1 v massive graph
       }
       else
       {
         double delta = 0;
         if (time_end > rna1.Teta) delta = rna1.Teta-time_start;
         else //time_end > Teta
           delta = time_end-time_start;
        // grafik spiralei
        for (int j = 0; j < exist_helises.size(); j++)
        {
          Integer I = (Integer) exist_helises.get(j);
          helix helix1 = (helix) this.get(I.intValue());
          helix1.graph[rna1.graph_step-1] += delta;
        }
        break;
        }
        if (next_time > time_end) // esli konec zhizni lezhit v predydushem elemente
        {
          // sootvetstvuyuschii otrezok vremeni dobavlyaem k elementu i
          // k grafiku spiralei
          for (int j = 0; j < exist_helises.size(); j++)
          {
            Integer I = (Integer) exist_helises.get(j);
            helix helix1 = (helix) this.get(I.intValue());
            helix1.graph[i] = helix1.graph[i] + (time_end - time_start);
          }
        }
        else // esli konec zhizni lezhit dal'she
        {
          // sootvetstvuyuschii otrezok vremeni dobavlyaem k elementu i
          // v grafike spiralei
          for (int j = 0; j < exist_helises.size(); j++)
          {
            Integer I = (Integer) exist_helises.get(j);
            helix helix1 = (helix) this.get(I.intValue());
            helix1.graph[i] = helix1.graph[i] + (next_time - time_start);
          }
        // i sdvigaem vremya nachala
        time_start = next_time;
        }
      }
    }


/*  //===================================
  // uporyadochivaem vector helices po nachalu spirali ls1
  void sort()
  {
    for (int i = 0; i < this.size(); i++)
    {
      helix sh = (helix) this.get(i);
      for (int j = i + 1; j < this.size(); j++)
      {
        helix sh1 = (helix)this.get(j);
        if (sh.ls1 > sh1.ls1)
        {
          this.remove(j);
          this.add(i, sh1);
        }
      }
    }
    int t = -1;
    for (int i = 0; i < this.size(); i++)
    {
      helix sh = (helix)this.get(i);
      int ls = sh.ls1;
      for (int j = i + 1; j < this.size(); j++)
      {
        helix sh1 = (helix)this.get(j);
        if (sh1.ls1 < ls)
        {
          ls = sh1.ls1;
          t = j;
        }
      }
      if (t != -1)
      {
        helix sh1 = (helix)this.get(t);
        this.remove(t);
        this.add(i, sh1);
        t = -1;
      }
      else continue;
    }
  }*/
}
