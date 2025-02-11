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


//=================================================

public class filter
{
    final int minhelix = 3; // razmer minimalnoi spirali

    //added by A.Favorov
    final int minpartialhelix = 3; // razmer minimalnogo kuska spirali

    // torchashego iz skolziashei spirali
    final int minloop = 3; // razmer minimal'noi petli

    //added by A.Favorov
    final double min_helix_living_time = 0.00002;

    // razreshaet/zapreschaet koncevye GU-pary
    final boolean if_GU = true;

    int Compatible[][]; // massiv sovmestimosti spiralei

    // nesovmestimy - (-1), sovmestimy - 0,
    // chastichno sovmestimy - nomer spirali so skol'z. petlei v vektore slid_heliset
    int heliset_size; // razmer vectora s normal'nymi spiralyami

    Vector heliset = new Vector(); // vse vozmozhnye spirali
    Vector from_sec_struct = new Vector(); // spirali iz struktury
    Vector helix_index = new Vector(); // nomera spiralei iz structury
    boolean[] exists = null; //
    boolean if_correct_structure = false; // govorit pravil'nuyu structuru narisoval Rfam ili net
    String curr_seq, curr_struct; // rassmatrivaemye posledovatel'nost' i struktura

    //========================================================
    // Kostruktor
    public filter(String seq, String structure)
    {
      curr_seq = seq;
      curr_struct = structure;
    //  find_helises(); // nahodim vse vozmozhnye spirali
      read_sec_struc(); // vytaskivaem spirali iz vtorichnoi struktury
    }

    //==============================================
    // vytaskivaem spirali iz skobochno-tochechnoi zapisi structury
    void read_sec_struc()
    {
      String struct = curr_struct; // rabochaya structura, iz kotoroi budem udalyat' naidennye spirali
      boolean correct_structure = true; // pokazyvaet vnutri plecha my ili net
      exists = new boolean[rna1.heliset1.size()];
      for (int i = 0; i < rna1.heliset1.size(); i++)
      {
        exists[i] = false;
      }
      int open_count = 0, close_count = 0;
      for (int i = 0; i < struct.length(); i++)
      {
        char ch = struct.charAt(i);
        if (ch =='(')
         {
           open_count++;
           continue;
         }
         if (ch == ')')
         {
           close_count++;
           continue;
         }
      }
      if (open_count!=close_count)
       {
         correct_structure = false;
       //  System.err.println("Not correct stucture");
        // return;
       }
      for (int i = 0; i < struct.length(); i++)
      {
        char ch = struct.charAt(i);
        int A = -1;
        int B = -1;
        int C = -1;
        int D = -1;
        if (ch == ')')
        {
          int j = -1;
          for (int j1 = i-1; j1 >= 0; j1--)
          {
            char ch_j = struct.charAt(j1);
            if (ch_j == '(')
            {
              j = j1;
              break;
            }
          }
          B = j;
          C = i;
          if (B != -1 && C != -1)
          for (int k = 0; k < struct.length(); k++)
          {
            if ((struct.charAt(j - k) == '(' && struct.charAt(i + k) == ')') &&
                ((j - k )!= 0 && (i + k) != struct.length()-1) ) continue;
            else
            {
              A = j - k + 1;
              D = i + k - 1;
              if ((j - k) == 0 || (i + k) == struct.length()-1 && (struct.charAt(j - k) == '(' && struct.charAt(i + k) == ')'))
              {
                A = j - k;
                D = i + k;
              }
              if ( A !=-1 && D != -1)
              {
                helix helix2 = new helix(A, B, C, D, curr_seq);
                from_sec_struct.add(helix2);
                // ischem nomer etoi spirali v heliset
                int n = rna1.heliset1.contain(helix2);
                if (n!=-1) exists[n] = true;
                else
                  {
                    System.err.println("Helix does not find");
                /*    System.err.println(rna1.ID + "\n" + "A = " + A +", B = " + B +
                                       ", C = " + C + ", D = " + D);//*/
                  }
                // pishem tochki vmesto naidennoi spirali
                String replace_shoulder = "";
                for (int m = 0; m < k; m++) replace_shoulder = replace_shoulder.
                    concat("."); // stroka iz tochek po dline spirali
                struct = struct.substring(0, A).concat(replace_shoulder).concat(struct.substring(B + 1, C)).
                    concat(replace_shoulder).concat(struct.substring(D + 1));
                break;
              }
            }
          }
        }
      }
      if (struct.length() != curr_struct.length())
       {
 //        System.err.println("Different size");
       }
      if (!correct_structure) // neodinakovoe kolichestvo "(" i ")"
      {
        String corr_struct = "";
        // delaem stroku iz tochek
        for (int i = 0; i < curr_seq.length(); i++)
          corr_struct = corr_struct.concat(".");
        // vstavlyaem tuda pravil'nye spirali
        for (int i = 0; i < from_sec_struct.size(); i++)
        {
          helix helix1 = (helix) from_sec_struct.get(i);
          String left = "", rigth = "";
          for (int j = 0; j < helix1.Left.length(); j++)
          {
            left = left.concat("(");
            rigth = rigth.concat(")");
          }
          corr_struct = corr_struct.substring(0,helix1.A).concat(left).concat(corr_struct.substring(helix1.B + 1, helix1.C)).
                  concat(rigth).concat(corr_struct.substring(helix1.D + 1));
        }
        curr_struct = corr_struct;
        //correct_structure = true;
        if_correct_structure = true;
      }
    }

    //=========================================
  // zapis' struktury v indeksah sparivaniya
  // esli ne sparennaya poziciya, to -1
  // esli sparennaya, to s kakoi
  int[] pair_index()
  {
    int sequence_length = curr_seq.length();
    int[] output_folding = new int[sequence_length]; // massiv
    for (int i = 0; i < sequence_length; i++) // vezde propisyvaem -1
    {
      output_folding[i] = -1;
    }
    for (int i = 0; i < from_sec_struct.size(); i++)
    {
      helix helix1 = (helix) from_sec_struct.get(i); // vytaskivaem sootvetstvuyuschuyu spiral'
      int k = 0;
      for (int j = helix1.A; j <= helix1.B; j++)
      { // levoe plecho
        output_folding[j] = helix1.D - k;
        k++;
      }
      k = 0;
      for (int j = helix1.C; j <= helix1.D; j++)
      { // pravoe plecho
        output_folding[j] = helix1.B - k;
        k++;
      }
    }
    return output_folding;
  }



    //============================================================
  // proveryaem sootvetstvie i korrektiruem vtor. strukturu
  String check()
  {
    String struct = curr_struct; // budem korrektirovat', esli nado
    // proveryaem sootvetstvie
    helix[] heliset_array = (helix[]) heliset.toArray(new helix[0]);
    helix[] sec_struct_array = (helix[]) from_sec_struct.toArray(new helix[0]);
    boolean correct_structure = true;
    for (int j = 0; j < sec_struct_array.length; j++)
    {
      if (correct_structure == true)
      {
        correct_structure = false;
        helix helix_j = sec_struct_array[j];
        int A2 = helix_j.A;
        int B2 = helix_j.B;
        int C2 = helix_j.C;
        int D2 = helix_j.D;
        for (int i = 0; i < heliset_array.length; i++)
        {
          int A1 = heliset_array[i].A;
          int B1 = heliset_array[i].B;
          int C1 = heliset_array[i].C;
          int D1 = heliset_array[i].D;
          if (A1 == A2 && B1 == B2 && C1 == C2 && D1 == D2)
          {
            correct_structure = true;
            break;
          }
        }
        // ne nashli podhodyaschei spirali, nado korrektirovat'
        //===============================
        // esli dlina spirali < minhelix, to udalyaem
        if (helix_j.Left.length() < minhelix)
        {
          String replace_shoulder = "";
          for (int m = 0; m < helix_j.Left.length(); m++) replace_shoulder =
              replace_shoulder.concat("."); // stroka iz tochek po dline spirali
          // udalyaem etu nepravil'nuyu spiral'
          struct = struct.substring(0, A2).concat(replace_shoulder).concat(struct.substring(B2 + 1, C2)).
              concat(replace_shoulder).concat(struct.substring(D2 + 1));
          correct_structure = true;
          if_correct_structure = true;
          continue;
        }

        for (int i = 0; i < heliset_array.length; i++)
        {
          helix helix_i = heliset_array[i];
          int A1 = helix_i.A;
          int B1 = helix_i.B;
          int C1 = helix_i.C;
          int D1 = helix_i.D;
          //===================================
          // esli helix_i yavlyaetsya podspiral'yu (prvil'naya spiral' yavlyaetsya podspiral'yu)
          if (A1 >= A2 && B1 <= B2 && C1 >= C2 && D1 <= D2)
          {
            String left = "", rigth = "", replace_shoulder = "";
            for (int k = 0; k < helix_i.Left.length(); k++)
            {
              left = left.concat("(");
              rigth = rigth.concat(")");
            }
            for (int m = 0; m < helix_j.Left.length(); m++) replace_shoulder = replace_shoulder.concat("."); // stroka iz tochek po dline spirali
            // udalyaem etu nepravil'nuyu spiral'
            struct = struct.substring(0, A2).concat(replace_shoulder).concat(struct.substring(B2 + 1, C2)).
                  concat(replace_shoulder).concat(struct.substring(D2 + 1));
            // i dobavlyaem pravil'nuyu
            struct = struct.substring(0, A1).concat(left).concat(struct.substring(B1 + 1, C1)).
                  concat(rigth).concat(struct.substring(D1 + 1));
            correct_structure = true;
            if_correct_structure = true;
            break;
          }
          //===================================
          // esli helix_j yavlyaetsya podspiral'yu helix_i (ostavlyaem kak est')
          if (A1 <= A2 && B1 >= B2 && C1 <= C2 && D1 >= D2)
          {
            correct_structure = true;
            if_correct_structure = true;
            break;
          }
        }
      }
      else break;
    }
    return struct;
  }

  //=====================================
    // nahodit vse vozmozhnye spirali
    void find_helises()
    {
      String seq = curr_seq; // tekushaya posledovatelnost
      int sl = seq.length();
      int[][] match = new int[sl][sl]; // array of matches for letters
      for (int i = 0; i < sl; i++)
      {
        for (int j = sl - 1; j > i; j--)
        {
          char chi = seq.charAt(i), chj = seq.charAt(j);
          // proveruaem yavlyaetsya li para (chi,chj) kanonicheskoi ili net
          if ( ( (chi == 'a') && (chj == 'u')) || ( (chi == 'u') && (chj == 'a')) ||
              ( (chi == 'c') && (chj == 'g')) || ( (chi == 'g') && (chj == 'c')))
          {
            match[i][j] = match[j][i] = 2;
          }
          if ( ( (chi == 'u') && (chj == 'g')) || ( (chi == 'g') && (chj == 'u')))
          {
            if (if_GU) match[i][j] = match[j][i] = 2; // razreshaem koncevye GU
            else match[i][j] = match[j][i] = 1; // zapreschaem
            // System.out.print(match[i][j]);
          }
        }
        //System.out.print("\n");
      }
  /*    // pechataem matricu match
      for (int i = 0; i < sl; i++)
      {
        System.out.print(s.charAt(i));
        for (int j = 0; j < sl; j++)
        {
          System.out.print(match[i][j]);
        }
        System.out.print("\n");
      } //*/

      // ischem obratnye diagonali v match
      int ls1 = -1, ls2 = -1, rs1 = -1, rs2 = -1;
      for (int i = (2 * minhelix - 1); i < sl; i++) // i - stolbec, j - stroka
      {
        int counter = 0; // schitaet 1 i 2
        for (int j = 0; 2 * j <= i; j++) // idem po obratnoi diagonali do osnovnoi
        {
          if ( ( (match[i - j][j] == 0) && (counter == 0)) || // net nachala dlya spirali
              ( (match[i - j][j] == 1) && (counter == 0)))continue; // ili nachinaetsya s GU
          if ( (match[i - j][j] == 2) && (counter == 0))
          {
            ls1 = j; // nachalo levogo plecha
            rs2 = i - j; // konec pravogo plecha
            counter++;
            continue;
          }
          else
          if ( ( (match[i - j][j] == 2) || (match[i - j][j] == 1)) &&
              (counter != 0))
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
                heliset.add(new helix(ls1, ls2, rs1, rs2, seq));
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
                      heliset.add(new helix(ls1, ls2, rs1, rs2, seq));
                      ls1 = ls2 = rs1 = rs2 = -1;
                      counter = 0;
                    }
                    if (match[i - j + t][j - t] == 1)
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
      // rassmatrivaem kusok matricy pod obratnoi diagonal'yu

      for (int j = 1; j < sl; j++)
      {
        int i = sl - 1;
        int counter = 0;
        for (int k = 0; k < i - j; k++)
        {
          if ( ( (match[i - k][j + k] == 0) && (counter == 0)) || // net nachala dlya spirali
              ( (match[i - k][j + k] == 1) && (counter == 0)))continue; // ili nachinaetsya s GU
          if ( (match[i - k][j + k] == 2) && (counter == 0))
          {
            ls1 = j + k; // nachalo levogo plecha
            rs2 = i - k; // konec pravogo plecha
            counter++;
            continue;
          }
          else
          if ( ( (match[i - k][j + k] == 2) || (match[i - k][j + k] == 1)) &&
              (counter != 0))
          {
            counter++;
            continue;
          }
          if ( (match[i - k][j + k] == 0) && (counter != 0)) // t.e. match[i][j]==0
          {
            if (counter >= minhelix)
            {
              if (match[i - k + 1][j + k - 1] == 2) // i zakanchivaetsya na GC ili AU
              {
                ls2 = j + k - 1; // konec levogo plecha
                rs1 = i - k + 1; // nachalo pravogo plecha
                heliset.add(new helix(ls1, ls2, rs1, rs2, seq));
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
                      heliset.add(new helix(ls1, ls2, rs1, rs2, seq));
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
      for (int i = 0; i < heliset.size(); i++)
      {
        helix sh = (helix)heliset.get(i);
        if ( (sh.A <= sh.C) && (sh.C <= sh.B))
        {
          heliset.remove(i--);
        }
      }

      // udalyaem odinakovye shpilki
      for (int i = 0; i < heliset.size(); i++)
      {
        helix sh = (helix)heliset.get(i);
        for (int j = i + 1; j < heliset.size(); j++)
        {
          helix sh1 = (helix)heliset.get(j);
          if ( (sh.A == sh1.C) && (sh.B == sh1.D) && (sh.C == sh1.A) &&
              (sh.D == sh1.B))
          {
            heliset.remove(j--);
          }
        }
      }

      // shpilki bez petli ili s petlei < 3bp, vida 3-5:8-6
      // nado smotret' mozhno li rasplesti nemnogo spiral', chtoby sdelat' normal'nuyu shpil'ku
      // esli nel'zya, to udalyaem
      for (int i = 0; i < heliset.size(); i++)
      {
        helix sh = (helix)heliset.get(i);
        if (sh.C - sh.B <= minloop)
        {
          int helix_length = sh.Left.length();
          if (helix_length == minhelix) // esli dlina spirali = minimal'noi dline, to nel'zya rasplesti
            heliset.remove(i--);
          else // esli dlina spirali bol'she minimal'noi i mozhno poprobovat' rasplesti
            for (int j = 1; j < helix_length; j++)
            {
              helix_length--;
              if (helix_length >= minhelix) // t.e. mozhno rasplesti
              {
                if (match[sh.B - j][sh.C + j] == 2) // esli teper' spiral' zakanchivaetsya na kanonicheskuyu paru
                {
                  sh.B = sh.B - j;
                  sh.C = sh.C + j;
                  sh.Left = seq.substring(sh.A, sh.B + 1); // levoe plecho
                  sh.Right = sh.reverse(seq.substring(sh.C, sh.D + 1)); // perevernutoe pravoe
                  sh.energy = sh.energy_calculation(sh.Left, sh.Right);
                  break;
                }
                else continue; // esli zakanchivaetsya na (g,u)
              }
              else heliset.remove(i--);
            }
        }
      }
      sort(); // uporyadochivaem vector helices

    //  return heliset;
      // udalyaem malozhivushe shpilki
      //added by A.Favorov
/*      for (int i = 0; i < this.size(); i++)
      {
        helix sh = (helix)this.get(i);
        if (sh.k_break > (1. / min_helix_living_time))
          this.remove(i--);
      }
      //System.out.println("");
      sort(); // uporyadochivaem vector helices
      heliset_size = this.size();

      //=====================================
      // ischem spirali, kotorye mogut obrazovat' structuru so skol'zhyaschei petlei
      // i kladem ih v slid_heliset
  /*    slid_heliset = new Vector();
      slid_heliset.add(0, null); // kladem pustoi element, chtoby indeksy u sk. spiralei byli >0
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
          if ( (A2 - A1 < D1 - D2) && // uslovia naklona vlevo
              (D1 >= C2) && (C2 > C1) && // pravye plechi peresekayutsya
              //           (B2 < A1) && // a levye ne peresekayutsya
              (C2 - C1 >= minpartialhelix) && // bez peresecheniya ostaetsya sparennymi ne men'she (minhelix-1) par
              (D2 - D1 >= minpartialhelix))
          {
            slid_heliset.add(new slid_helix(i, j, helix1, helix2));
            continue;
          }
          if ( (A2 - A1 > D1 - D2) && // naklon vpravo
              (B1 > B2) && (B2 >= A1) && // levye plechi peresekayutsya
              //           (C2 > D1) && // a pravye net
              (B1 - B2 >= minpartialhelix) &&
              (A1 - A2 >= minpartialhelix))
          {
            slid_heliset.add(new slid_helix(i, j, helix1, helix2));
            continue;
          }
        }
      }*/

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
        }//
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
        for (int i = 1; i < slid_heliset.size(); i++)
        {
          slid_helix sh = (slid_helix) slid_heliset.get(i);
          Compatible[sh.slip_parts[0]][sh.slip_parts[1]] =
              Compatible[sh.slip_parts[1]][sh.slip_parts[0]] = i;
        }*/
    }

//===============================================
    // proverka shpilek na sovmestimost
    public int isCompatible(helix a, helix b)
    {
      int result = 0;
      if ( (a.B <= b.A) && (b.A <= a.D) && (a.C <= b.D)) result = -1; // psevdouzel
      else
      if ( (b.B <= a.A) && (a.A <= b.C) && (b.C <= a.D)) result = -1;
      else
      if ( ( (a.A <= b.A) && (b.A <= a.B)) || ( (b.A <= a.A) && (a.A <= b.B)))
        result = -1;
      else
      if ( ( (b.C <= a.D) && (a.D <= b.D)) || ( (a.C <= b.D) && (b.D <= a.D)))
        result = -1;
      return result;
    }

//===================================
    // uporyadochivaem vector helices po koncu spirali rs2
    void sort()
    {
      for (int i = 0; i < heliset.size(); i++)
      {
        helix sh = (helix)heliset.get(i);
        for (int j = i + 1; j < heliset.size(); j++)
        {
          helix sh1 = (helix)heliset.get(j);
          if (sh.D > sh1.D)
          {
            heliset.remove(j);
            heliset.add(i, sh1);
          }
        }
      }
      int t = -1;
      for (int i = 0; i < heliset.size(); i++)
      {
        helix sh = (helix)heliset.get(i);
        int ls = sh.D;
        for (int j = i + 1; j < heliset.size(); j++)
        {
          helix sh1 = (helix)heliset.get(j);
          if (sh1.D < ls)
          {
            ls = sh1.D;
            t = j;
          }
        }
        if (t != -1)
        {
          helix sh1 = (helix)heliset.get(t);
          heliset.remove(t);
          heliset.add(i, sh1);
          t = -1;
        }
        else continue;
      }
    }

    //===================================================
    // zapolnyaem grafik dlya spiralei
/*    void add_to_graph(double CurrTime, double LiveTime, Vector exist_helises)
    {
      // ischem mesto v massive graph
      double time_start = CurrTime; //nachalo zhizni struktury
      double time_end = CurrTime + LiveTime;
      if (CurrTime + LiveTime > rna1.Teta) time_end = rna1.Teta; // konec zhizni
      int n_start = (int) ( (rna1.graph_step * time_start) / rna1.Teta); // element massiva graph, kuda popadaet nachalo
      int n_end = (int) ( (rna1.graph_step * time_end) / rna1.Teta); // element massiva graph, kuda popadaet konec
      for (int i = n_start; i <= n_end; i++)
      {
        double next_time = 0;
        if (i + 1 < rna1.graph_step)
        {
          next_time = ( (i + 1) * rna1.Teta) / rna1.graph_step; // vremya, sootvetstvuyuschee elementu i+1 v massive graph
        }
        else
        {
          double delta = 0;
          if (time_end > rna1.Teta) delta = rna1.Teta - time_start;
          else //time_end > Teta
            delta = time_end - time_start;
            // grafik spiralei
          for (int j = 0; j < exist_helises.size(); j++)
          {
            Integer I = (Integer) exist_helises.get(j);
            helix helix1 = (helix)this.get(I.intValue());
            helix1.graph[rna1.graph_step - 1] += delta;
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
            helix helix1 = (helix)this.get(I.intValue());
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
            helix helix1 = (helix)this.get(I.intValue());
            helix1.graph[i] = helix1.graph[i] + (next_time - time_start);
          }
          // i sdvigaem vremya nachala
          time_start = next_time;
        }
      }
    }*/

    //===================================
 /*   // uporyadochivaem vector helices po nachalu spirali ls1
    void sort()
    {
      for (int i = 0; i < this.size(); i++)
      {
        helix sh = (helix)this.get(i);
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

