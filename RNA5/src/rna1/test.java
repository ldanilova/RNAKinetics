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


// schitaet vse sovpadeniya i nesovpadeniya mezhdu paroy struktur, i ocenki QL, QQ i CC
class estimation
  {
    int PC = 0, PW = 0, PU = 0, UC = 0, UW = 0, L = 0;
    double QL = 0., QQ = 0., CC = 0.;

    public estimation(int[] index_foung_struct, int[] index_pattern)
    {
      L = index_foung_struct.length;
      for (int j = 0; j < L; j++)
      {
        int predict = index_foung_struct[j];
        int pattern = index_pattern[j];
        //PC (paired correctly)
        if (predict == pattern && predict > -1)
        {
          PC++;
          continue;
        }
        // PW (paired wrongly)
        if (predict > -1 && pattern > -1 && predict != pattern)
        {
          PW++;
          continue;
        }
        // PU (paired to unpaired)
        if (predict > -1 && pattern == -1)
        {
          PU++;
          continue;
        }
        // UC (unpaired correctly)
        if (predict == pattern && predict == -1)
        {
          UC++;
          continue;
        }
        // UW (unpaired wrongly)
        if (pattern > -1 && predict == -1)
        {
          UW++;
          continue;
        }
      }
      double length = L;
      QL = (PC+UC+PU)/length;
      QQ = PC/(PW+PU+PC+UW);
      CC = (PC *(PU+UC) - PW * UW) / Math.sqrt( (PC + PW)*( (PU+UC)+ UW)*( PC + UW)*( (PU+UC)+ PW) );
    }
  }




public class test
{


  public static void main(String[] args)
  {
    ReadFile readFile1 = new ReadFile("strong_helices1.txt"); // file with input dataset
    String[] sequences = readFile1.get_sequences();
    String[] names = readFile1.get_names();
    String[] structures = readFile1.get_structures();
    int[] index_foung_struct; // moya naidennaya strucktura
    WriteFile writelog = new WriteFile("logfile.log", false);
    WriteFile outputfile = new WriteFile("result.txt", false);
    for (int i = 0; i < names.length; i++)
    {
      long t1 = System.currentTimeMillis();
      foldings st1 = null;

      rna1 rna = new rna1(sequences[i], names[i]);
      int length = sequences[i].length();
      int result = rna.elongate();
      if (result == -1)// esli programma rabotaet bolshe nekotorogo vremeni
      {
        outputfile.Writeln(names[i] + "\t" + "works too long");
        continue;
      }
      if (result == 0) // esli net spiralei
      {
        outputfile.Writeln(names[i] + "\t" + "There is no helices");
        continue;
      }
      if (result == 1) // esli spiral' tol'ko odna
      {
        boolean exists[] = new boolean[1];
         exists[0] = true;
         foldings folding = new foldings(exists, rna.heliset1);
         index_foung_struct = folding.pair_index(length);
        // outputfile.Writeln(names[i] + "\t" + "There is only one helix");
       }
      else // esli vse v poryadke i spiralei bol'she 1
      {
        // vytaskivaem element s luchshei naidennoi strukturoi
        st1 = rna.get_best_structure();
        // vytaskivaem iz nee skobochni tochechnuyu zapis'
        //String found_structure = st1.folding.bracket_dots(sequences[i].length());
        // i indeksnuyu zapis'
        index_foung_struct = st1.pair_index(length); // moya strucktura
        // vytaskivaem iz otveta spirali
      }
      //========================================
      // sravnivaem otvety
      // ocenochnyt parametry dlya struktury, pervoi po spisku
      int PC = 0, PW = 0, PU = 0, UC = 0, UW = 0;
      estimation e = null;
      // ocenochnye parametry dlya pravil'noi struktury
     int PC1 = 0, PW1 = 0, PU1 = 0, UC1 = 0, UW1 = 0;
      estimation e1 = null;
      String comment = "";

      filter filter1 = new filter(sequences[i], structures[i]);
      // zapisyvaem v indexnoi zapisi
      int[] index_pattern = filter1.pair_index(); // Rfam
      // schitaem vse sovpadeniya i nesovpadeniya, i ocenki QL, QQ i CC
      e = new estimation(index_foung_struct, index_pattern);

      Comparator_exists com = new Comparator_exists();
      if (com.compare(filter1.exists, st1.exists) != 0)
      {// esli pervay struktura nepravilnaya
        // proveryaem, mozhet byt' ona nizhe po spisku
        foldings st2 = rna.contain_structure(filter1.exists);
        if (st2 != null)
        {
          int[] index_struct = st2.pair_index(length);
          e1 = new estimation(index_struct, index_pattern);
          comment = "Below in the list";
        }
        else
        // esli ne nashli v spiske, izhem v group_list

        if (rna.groups_list1.hash.containsKey(filter1.exists))
        {// esli est' v spiske grupp
          // esli est', vytaskivaem nomer gruppy
          Integer group_index = (Integer) rna.groups_list1.hash.get(
              filter1.exists);
          // i samu gruppu po etomu nomeru
          group gr = (group) rna.groups_list1.get(group_index.intValue());
          if (!gr.completed)
            {
              System.err.println("Group is not completed");
              comment = "Group is not completed. ";
            }
          else
          {
            int n = gr.find_element(filter1.exists);
            if (n == -1) System.err.println("Element does not find");
            else
            {
              int[] index_struct = gr.elements[n].pair_index(length);
              e1 = new estimation(index_struct, index_pattern);

              comment = comment.concat("In group #" + group_index.intValue() + " element #" + n);
              if (gr.find_element(st1.exists) != -1) comment = comment.concat(". The first group in the list");
            }
          }
        }
        else // net v group liste
        comment = "Not in the group list";
      }
      PC = e.PC;
      PW = e.PW;
      PU = e.PU;
      UC = e.UC;
      UW = e.UW;

      //==============================
      // zapisyvaem otvet
      if (result == 1) // esli spiral' tol'ko odna
      {
        outputfile.Writeln(names[i] + "\t" + PC + "\t" + PW + "\t" +
                         PU + "\t" + UC + "\t" + UW + "\t" +
                         length + "\t" + /* "\t" +
                         PW_corr + "\t" + UW_corr + "\t" + seq_length_corr + "\t" +*/
                         "There is only one helix");
      }
      else
        if (e1 == null)
      outputfile.Writeln(names[i] + "\t" + PC + "\t" + PW + "\t" +
                         PU + "\t" + UC + "\t" + UW + "\t" +
                         length + "\t" + comment/*"\t" +
                         PW_corr + "\t" + UW_corr + "\t" + seq_length_corr*/);
    else
    {
      PC1 = e1.PC;
      PW1 = e1.PW;
      PU1 = e1.PU;
      UC1 = e1.UC;
      UW1 = e1.UW;
      outputfile.Writeln(names[i] + "\t" + PC + "\t" + PW + "\t" +
                               PU + "\t" + UC + "\t" + UW + "\t" +
                           length + "\t" + "\t" +
                              PC1 + "\t" + PW1 + "\t" + PU1 + "\t" + UC1 + "\t" + UW1 + "\t" +
                               "\t" + comment );
    }

      long t2 = System.currentTimeMillis();
      writelog.Writeln(names[i] + "\t" + (t2 - t1));
    }
  }

}
