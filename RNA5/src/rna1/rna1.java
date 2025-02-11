package rna1;

import java.util.*;
//import java.io.*;

/**
 * <p>Title: ReadFile</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Luda
 * @version 1.0
 */

//$Id: rna1.java,v 1.7 2003/10/29 01:35:55 favorov Exp $
//======================================================
// Svorachivaem RNK po mere rosta cepi i hranim vse promezhutochnye sostoyaniya
// rassmatrivaem nesovmestimost'
// Konsol'nyi variant
//======================================================

public class rna1
{
  static double k_gr = 40;// konstanta skorosti rosta cepi, t.e. poyavleniya ocherednogo nukleotida (c^-1)
  static double k0 = 1.e+7; // konstanta skorosti obrazovaniya ocherednoi pary (c^-1)
  static double Teta = 20;// process prekraschaetsya, kogda time stanet ravno Teta (c)
  static int M = 1; // kolichestvo prohodov (run)
  static double threshold_dis = 1.e+4; // porog na kostantu pri ob'edinenii v foldigs v gruppy
  static double threshold_form = 1.e+4; // porog na kostantu pri ob'edinenii v foldigs v gruppy
  static double group_perehod = 1.e-4; // porog na konstantu pri perehode mezhdu gruppami
  static String ID = "test";// id zadachi, peredavaemoe iz web-form
  static int graph_step = 100; // razmer shaga vremeni dlya grafikov 1 i 2 (Teta/graph_step)
  static int graph_seq = 5; // kolichestvo srezov vremeni dlya grafika 3
  static Random Random1 = new Random();

  //static int count_lifeperiod = 0; // schitaet kolichestvo new LifePeriod
  static double grown_time = 0; // moment vremeni kogda vyrosla vsya posledovatel'nost'
  static boolean print_graph = false; // nado li pechatat' graphiki
  static boolean print_only_helix_index = true; // nado li pechatat' spirali v polnoi forme v strukturah
  static boolean for_site = true; // kompeliruetsya variant programmy dlya saita ili net

  // parametry dlya var'irovaniya energii
  double helix_energy_var = 0; // procent var'irovaniya dlya energii stekinga (0 - na 0%; 0.2 - na +-10%; 0.4 - na +-20%)
  double loops_energy_var = 0; // procent var'irovaniya dlya energii petel'
  static double[] energy_helix_coef = new double[36]; // koefficienty dlya varirovaniya energii stekinga
  static double[] energy_loops_coef = new double[37]; // koefficienty dlya varirovaniya energii petel'

  // parametry dlya heliset
  static int minhelix = 3; // razmer minimalnoi spirali
  static int minpartialhelix = 3; // razmer minimalnogo kuska spirali
  static double helix_threshold = -7000; // porog otsecheniya spiralei
  static boolean if_GU = true; // razreshaet/zapreschaet koncevye GU-pary

  static int min_stack = 2; // min kolichestvo sparennyh nukleotidov v strukture so skol'zyaschei petlei

  static WriteFile writestruct; // file with perehodami

//  boolean exception = false; // esli programma vydast soobschenie o nehvatke pamyati
  static String sequence = "";
  String outputfile = "";
  TreeMap st_time; // sostonyaniya sortirovannye po vremeni
  //TreeMap states; // sostoyaniya s perehodami, konstantami perehoda i schetchikom
  Comparator_exists C_exists; // klyuch, po kotoromu budut uporyadocheny elementy states
  Comparator_time C_time; // klyuch' dlya sortirovki po vremeni i energii
  TreeMap st_energy; // sostonyaniya sortirovannye po energii
  Comparator_energy C_energy;
  // vyhodnoi file
  WriteFile writeoutput; // file with output data
  static heliset heliset1;
  static groups_list groups_list1; // gruppy sostoyanii
  Comparator_vector C_vector; // klyuch dlya loops_list
  static TreeMap loops_list; // spisok petel'

  //static WriteFile writestruct; // file with perehodami


  //============================================
  // KONSTRUKTOR

  public rna1(String seq, String name)
  {
    sequence = seq;
    if (!for_site) // esli ne dlya saita
    ID = name.replace('/','_');

    if (for_site) // esli dlya saita
    {
      print_graph = true; // nado pechatat' graphiki
      print_only_helix_index = true; // i tol'ko indeksy spiralei v strukture
//      if_GU = false;
    }

    //=========================================
    sequence = sequence.trim();
    sequence = sequence.toLowerCase();
    sequence = sequence.replace('t', 'u');
    // ReadFile readFile1 = new ReadFile("test1.txt");// file with input dataset
    if (outputfile != "")
    { // esli imya vyhocnogo file zadano
      writeoutput = new WriteFile(outputfile, false);
    }
    else { // esli ne zadano
      writeoutput = new WriteFile(ID + ".out", false); // file with output data
    }

    //==========================================
    // zapolnyaem massiv koeff. dlya varirovaniya energii stakinga i zapominaem ih v fail
//    WriteFile writecoeff = new WriteFile(ID + ".coef", false); // file with coefficients
//    writecoeff.Writeln("stacking energy coefficients");
    for (int i = 0; i < energy_helix_coef.length; i++)
    {
      energy_helix_coef[i] = (1 - helix_energy_var/2) + helix_energy_var * Random1.nextDouble();
//      writecoeff.Writeln(i + "\t"+ energy_helix_coef[i]);
    }
//   writecoeff.Writeln("loops energy coefficients");
    for (int i = 0; i < energy_loops_coef.length; i++)
    {
      energy_loops_coef[i] = (1-loops_energy_var/2) + loops_energy_var * Random1.nextDouble();
//      writecoeff.Writeln(i + "\t"+ energy_loops_coef[i]);
    }//*/
  }

  //================================================================
  // OBRABOTKA POSLEDOVATEL'NOSTI
  // resul'taty:
  // 0 - kogda net spiralei
  // 1 - kogda odna spiral'
  // 2 - kogda > 2
  // -1 - kogda schitaet slishkom dolgo

  int elongate()
  {
      //========================================
      // zapisyvaem parametry zadachi
      writeoutput.Writeln("ID = " + ID);
      writeoutput.Writeln(sequence);
//      writeoutput.Writeln("sequence length = " + sequence.length());
      writeoutput.Writeln("Time = " + Teta);
      writeoutput.Writeln("M = " + M);
      writeoutput.Writeln("k_growth = " + k_gr);
      writeoutput.Writeln("k_nucl = " + k0);
      writeoutput.Writeln("helix threshold = " + helix_threshold);
 //     writeoutput.Writeln("destroy threshold = " + threshold_dis);
  //    writeoutput.Writeln("formation threshold = " + threshold_form);
      writeoutput.Writeln("min helix length = " + minhelix);
      writeoutput.Writeln("allows final GU = " + if_GU);
      try
      {
        // nashli vse vozmozhnye spirali
        heliset1 = new heliset(sequence);
      }
      catch (Exception ex)
      {
        WriteFile writelog = new WriteFile(ID + ".log", false);
        writelog.Writeln("Heliset creation. Out of memory");
        System.exit(1);
      }

        int heliset_size = heliset1.heliset_size;
        // log file
        WriteFile writelog = new WriteFile(ID + ".log", false);
 //       writestruct = new WriteFile("perehody.out", false); // file with structures

        //================================
        // esli M = 0, to prosto raspechatyvaem spisok spiralei
        if (M == 0 && !heliset1.isEmpty())
        {
          //zapisyvaem vse spirali vozmozhnye v tekuzchei posledovatel'nosti
          writeoutput.Writeln("All helices");
          for (int i = 0; i < heliset_size; i++)
          {
            writeoutput.Writeln(i + "\n" + heliset1.get(i).toString());
          }
          writeoutput.Writeln("slip helices");
          for (int i = 1; i < heliset1.slid_heliset.length; i++)
          {
            slid_helix helix1 = (slid_helix) heliset1.slid_heliset[i];
            writeoutput.Writeln(i + "\n" + helix1.slip_parts[0] + ", " +
                                helix1.slip_parts[1] + "\n" +
                                helix1.toString() + "\n" + "A1=" + helix1.A1 +
                                " B1=" + helix1.B1 + " C1=" + helix1.C1 + " D1=" +
                                helix1.D1
                                + " A2=" + helix1.A2 + " B2=" + helix1.B2 +
                                " C2=" + helix1.C2 + " D2=" + helix1.D2 + "\n" +
                                "cross = " + helix1.cross + "\n" +
                                "cross energy = " + helix1.energy_cross + "\n" +
                                "energy1 = " + helix1.energy1 + "\n" +
                                "energy2 = " + helix1.energy2 + "\n");
          }

          return 2;
        }

        // smotrim est' li spirali
        switch (heliset_size)
        {
          case 0:
          { // esli net spiralei
            System.err.println("There is no helices");
            writelog.Writeln("There is no helices");
            writeoutput.Writeln("There is no helices");
            return 0;
          }
          case 1:
          {
            System.err.println("There is only one helix");
            writelog.Writeln("There is only one helix");
            writeoutput.Writeln("There is only one helix");
            boolean exists[] = new boolean[1];
            exists[0] = true;
            foldings folding = new foldings(exists, heliset1);
            writeoutput.Writeln("structure energy = " + folding.folding_energy);
            writeoutput.Writeln(folding.bracket_dots(sequence.length()));
            return 1;
          }
          default:
          { // kogda est' spirali

            C_exists = new Comparator_exists();
            C_time = new Comparator_time();
            C_energy = new Comparator_energy();
            st_time = new TreeMap(C_time);
            st_energy = new TreeMap(C_energy);

            groups_list1 = new groups_list();

            C_vector = new Comparator_vector();
            loops_list =  new TreeMap(C_vector);


            int m = 1;

            // vremya nachala structurnyh perestanovok
            long t1 = System.currentTimeMillis();

            structure structure1 = new structure(heliset1, sequence.length());

//         WriteFile writestates = new WriteFile("states.txt", false);

              //==========================================
              // cikl po kolichestvu prohodov M
              while (m <= M)
              {
                //========================================
                // zapisyvaem nomer tecuschego runa v log-file
                // log file
                writelog = new WriteFile(ID + ".log", false);
                writelog.Writeln("current run m = " + m); // zapisyvaem nomer prohoda
                writelog.Writeln("M = " + M); // zapisyvaem obschee chislo prohodov
                writelog.Writeln("start time = " + (long) t1 / 1000); // zapisyvaem tekuschee vremya
                writelog.Writeln("current time = " +
                                 (long) (System.currentTimeMillis() / 1000)); // zapisyvaem tekuschee vremya
         //*/
//                long t2 = System.currentTimeMillis();
                try
                {
                  structure1.run(groups_list1);
                }
                catch (OutOfMemoryError ex1)
                {
                  writelog.Writeln("Out of memory");
                  System.exit(1);
            //      break;
                }

 //               System.err.println("m = " + m + "\t" + (System.currentTimeMillis() - t2));
                // esli programma rabotaet bol'she chasa i eshce ne proshlo i poloviny, to nado vyhodit'
                //    if (t2 - t1 > 3600000 && m < M / 2)return -1;
                //     esli programma rabotaet bol'she 5 min, to nado vyhodit'
   //             if (!for_site) // esli ne dlya saita
   //               if (t2 - t1 > 300000)return -1;

                //writestates.Writeln(m + "\t" + states.size() + "\t" + structure1.count_new);

                m++;
              }

              writeoutput.Writeln("structure permutations # = " +
                                  structure1.GlobalCount);

              // zapisyvaem vremya, kogda vyrosla vsya posl-t', esli ona uspela vyrosti
              // ili poziciyu, esli ne uspela
 /*             if (grown_time != 0)
                writeoutput.Writeln("grown time = " + grown_time);
              else writeoutput.Writeln("grown position = " + structure1.CurPos);
   */         if (print_graph == true)
            {
              // zapisyvaem parametry dlya grafikov
 /*            writeoutput.Writeln("graph step for structures = " + graph_step);
              writeoutput.Writeln("graph step for sequence = " + graph_seq);
 */
             writeoutput.Writeln("step 1 = " + graph_step);
              writeoutput.Writeln("step 2 = " + graph_seq);// zapisyvaem grafiki dlya posledovatel'nosti
              for (int i = 0; i < graph_seq; i++)
              {
                StringBuffer sb = new StringBuffer();
                for (int j = 0; j < sequence.length(); j++)
                {
                  // delim na kolichestvo prohodov
                  structure1.graph_3[i][j] = (structure1.graph_3[i][j] / (M));
                  // delim na shag, chtoby poluchilas' veroyatnost'
                  structure1.graph_3[i][j] = (structure1.graph_3[i][j] /
                                              (Teta / graph_seq));
                  sb.append(structure1.graph_3[i][j] + ", ");
                }
                writeoutput.Writeln("sequence_graph " + i + ": " + sb.toString());
              }
            }
            long t2 = System.currentTimeMillis();
//            writeoutput.Writeln("work time = " + (t2 -t1));
            writeoutput.Writeln("");

            //================================================
            //zapisyvaem vse spirali vozmozhnye v tekuzchei posledovatel'nosti
            writeoutput.Writeln("All helices");
            for (int i = 0; i < heliset1.heliset_size; i++)
            {
              writeoutput.Writeln(i + "\n" + heliset1.get(i).toString());
            }
 /*           writeoutput.Writeln("slip helices");
            for (int i = 1; i < heliset1.slid_heliset.size(); i++)
            {
              slid_helix helix1 = (slid_helix) heliset1.slid_heliset.get(i);
              writeoutput.Writeln(i + "\n" + helix1.slip_parts[0] + ", " +
                                  helix1.slip_parts[1] + "\n" +
                                  helix1.toString() + "\n" + "A1=" + helix1.A1 +
                                  " B1=" + helix1.B1 + " C1=" + helix1.C1 + " D1=" +
                                  helix1.D1
                                  + " A2=" + helix1.A2 + " B2=" + helix1.B2 +
                                  " C2=" + helix1.C2 + " D2=" + helix1.D2 + "\n" +
                                  "cross = " + helix1.cross + "\n" +
                                  "cross energy = " + helix1.energy_cross + "\n" +
                                  "energy1 = " + helix1.energy1 + "\n" +
                                  "energy2 = " + helix1.energy2 + "\n");
            }//*/


            //=========================================
            // sortiruem group_list po vremeni i energii
            if (!groups_list1.isEmpty())
            {
              for (int i = 0; i < groups_list1.size(); i++) // i perebiraem ih iteraterom
              {
                group gr = (group) groups_list1.get(i); // izvlekaem sootvetstvuyuschii element
                if (gr != null && gr.completed)
                {
                  st_time.put(new Float(gr.live_time), gr);
                  foldings f = (foldings) gr.elements[gr.best];
                  st_energy.put(new Float(f.folding_energy), gr);
                }
              }

              //===================================
              // zapisyvaem statistiki uporyadochennye po vremeni
              writeoutput.Writeln("SORT BY TIME");
              // vsyu statistiku
              writeoutput.Writeln("All statistics");
              Set key_set = st_time.keySet(); // vytaskivaem set klyuchei
              for (Iterator iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
              {
                Float key = (Float) iter.next();
                group gr = (group) st_time.get(key);
                writeoutput.Writeln(gr.toString(sequence.length()));
              }

              WriteFile out = new WriteFile(ID + "_groups.out",false);
              // zapisyvaet tol'ko polnye gruppy
              out.Writeln(groups_list1.toString(true));
 //             writeoutput.Writeln("states.size = " + states.size());
              //==================================================
              // zapisyvaem statistiku uporyadochennye po energii
              /*         writeoutput.Writeln("SORT BY ENERGY");
                       //structure1.statistic1.write_statistics(writeoutput, structure1, false, true);
                       key_set = st_energy.keySet(); // vytaskivaem set klyuchei
                       for (Iterator iter = key_set.iterator(); iter.hasNext(); )  // i perebiraem ih iteraterom
                       {
                         Float key = (Float) iter.next();
                         state_element st1 = (state_element) st_energy.get(key);
                         writeoutput.Writeln(st1.toString(sequence.length()));
                       }//*/

               // zapisyvaem v log-file nomer poslednego prohoda, kak ukazanie, chto output-file zapisan
 /*             writelog = new WriteFile(ID + ".log", false);
              writelog.Writeln("current run m = " + M); // zapisyvaem nomer prohoda
              writelog.Writeln("M = " + M); // zapisyvaem obschee chislo prohodov
              writelog.Writeln("start time = " + (long) t1 / 1000); // zapisyvaem tekuschee vremya
              writelog.Writeln("current time = " +
                               (long) (System.currentTimeMillis() / 1000)); // zapisyvaem tekuschee vremya

             WriteFile writelist = new WriteFile("groups_list.out", false);
              int all_elements = 0; // schitaet skol'ko elementov vo vseh gruppah
              for (int i = 0; i < structure1.groups_list1.size(); i++)
              {
                //writelog.Writeln(i + "\n");
                group gr1 = (group) structure1.groups_list1.get(i);
                if (gr1 != null)
                {
                  all_elements += gr1.elements.length;
                  for (int j = 0; j < gr1.elements.length; j++)
                  {
                    StringBuffer sb = new StringBuffer();
                    sb.append(i + "\n" + gr1.elements[j].toString() + "\n");
                    // pechataem nomer po poryadku
                    for (int k = 0; k < heliset1.size(); k++)
                    {
                      sb.append(k + "\t");
                    }
                    sb.append("\n"); // pechataem exists
                    for (int k = 0; k < heliset1.size(); k++)
                    {
                      if (gr1.elements[j].base_folding.exists[k]) sb.append("t" +
                          "\t");
                      else sb.append("f" + "\t");
                    }
                    sb.append("\n");
                    // pechataem konstanty perehoda
                    for (int k = 0; k < gr1.elements[j].outside_const.length; k++)
                    {
                      sb.append( (int) gr1.elements[j].outside_const[k] + "\t");
                    }
                    sb.append("\n"); //
                    writelist.Writeln(sb.toString()); // zapisyvaem nomer prohoda
                  }
                }

              }
              writelist.Writeln("hash size = " + structure1.groups_list1.hash.size() + "\n" +
                               "all elements = " + all_elements);
              //*/


   /*                // zapisyvaem massiv so srtukturami i konstanty perehodov
                          WriteFile writestruct = new WriteFile(ID + "_structures.out", false);// file with structures
                          key_set = states.keySet(); // vytaskivaem set klyuchei
                          for (Iterator iter = key_set.iterator(); iter.hasNext(); )  // i perebiraem ih iteraterom
                          {
                            boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
                            state_element st1 = (state_element) states.get(key_exists); // izvlekaem sootvetstvuyuschii element
                            StringBuffer sb = new StringBuffer();
                            // pechataem nomer po poryadku
                            for (int i = 0; i< st1.perehod_amount.length; i++)
                            {
                              sb.append(i + "\t");
                            }
                            sb.append("\n");// pechataem exists
                            for (int i = 0; i< key_exists.length; i++)
                            {
                              if (key_exists[i]) sb.append("t" + "\t");
                                  else sb.append("f" + "\t");
                            }
                            sb.append("\n");
                            // pechataem konstanty perehoda
                            for (int i = 0; i< st1.konst.length; i++)
                            {
                              sb.append((int)st1.konst[i] + "\t");
                            }
                            sb.append("\n");

                            // pechataem kolichestvo perehodov
                            for (int i = 0; i< st1.perehod_amount.length; i++)
                            {
                              sb.append(st1.perehod_amount[i] + "\t");
                            }
                            sb.append("\n");
                            // pechataem summarnoe vremya, kotoroe strukrura byla v etom sostoyanii
                            sb.append("live time = " + st1.live_time/rna1.M);
                            writestruct.Writeln(sb.toString() + "\n");
                          }//*/
               // System.err.println("count new = " + structure1.count_new);
               // vremya konca structurnyh perestanovok
              t2 = System.currentTimeMillis();
  //             System.err.println("work time = " + (t2-t1));
            }
          }
        }

        return 2;
      }

      //=========================================================================
      // RAZBIRAET KOMANDNUYU STROKU

      public static void main(String[] args)
      {
        try
        {

          String sequence = "";
    //    double Teta, k_gr, k0, helix_threshold;
    //    int M, minhelix;
    //    boolean if_GU;
          rna1 rna = null;

          //==============================================
          // razbiraem komandnuyu stroku
          for (int i = 0; i < args.length; i++)
          {
            if (args[i].indexOf("-id") != -1)
            { // esli v etoi stroke soderzhitsya id zadachi
              i++;
              ID = args[i].trim(); // vytaskivaem id zadachi i udalyaem lishnie probely
              continue;
            }
            if (args[i].indexOf("--id") != -1)
            { // esli v etoi stroke soderzhitsya id zadachi
              i++;
              ID = args[i].trim(); // vytaskivaem id zadachi
              continue;
            }
            if (args[i].indexOf("-seq") != -1)
            { // esli posledovatel'nost'
              i++;
              sequence = args[i].trim();
              continue;
            }
            if (args[i].indexOf("--seq") != -1)
            { // esli posledovatel'nost'
              i++;
              sequence = args[i].trim();
              continue;
            }
            if (args[i].indexOf("-in") != -1)
            { // esli imya vhodnogo fila
              i++;
              String s = args[i].trim();
              ReadFile readFile1 = new ReadFile(s); // file with input dataset
              String[] sequences = readFile1.get_sequences();
              sequence = sequences[0];
              continue;
            }
            if (args[i].indexOf("--infile") != -1)
            { // esli imya vhodnogo fila
              i++;
              String s = args[i].trim();
              ReadFile readFile1 = new ReadFile(s); // file with input dataset
              String[] sequences = readFile1.get_sequences();
              sequence = sequences[0];
              continue;
            }
          }

//          sequence = "ggggcuauagcucagcugggagagcgccugcuuugcacgcaggaggucugcgguucgaucccgcauagcuccacca";
                    // (((((((..((((........)))).(((((.......)))))......((((.......)))).)))))))....	64=7:48, (1,4,11,15)

 sequence = "CGGCGTTGATAGGATTTAGTAGTATTTCGATTTCTGATTTCAGAGAGCTGGTGGTCGGTGCGAACCAGTACAGAGCGAATTATGAATTACCCCCTGGAGCTTCTTTTACGAAACGTAAGGAGTAGTGAAAGACGGTTATAGACCGTTATGTCTAAAGAGTGGTGAAACGAAACTGTTTCACAATTTAGGGTGGTACCGCGAATTTTTCGTCCCTGCATATATTGCAGGGGCGTTTTTATTTTATAAGCGCATATCATATGACCCTTTATTTTGTGATACATT"; //
//          sequence = "UAGUUACUGGGGGUGCCCGCUUUCGGGCUGAGAGAGAAGGCAAGCUUCUUAACCCUUUGGACCUGAUCUGGUUCGUACCAGCGUGGGGAAGUAGAGGAAUUGUUUUUGUUAUU";
//         sequence = "gccgggguggcccagccugguagggcgucggccugcuaagccgaugauccguuaaggaucgcgcggguucaaaucccguccccggcg"; //M32222.1_1277-136
//         sequence ="ccccctgcgcgaaaaacgcgctcacactttttgtgtgactctcaaaaagagagtggggg"; // iskusstvennaya tRNA

//          sequence ="cccccaaaaaagggggaaaaaagggggaaaaaaccccc";

//         ID = "X16886.1/923-85_corr"; // (((((((..((((......)))).(((((.......)))))....(((((.......)))))))))))).

          rna = new rna1(sequence, ID);

          for (int i = 0; i < args.length; i++)
          {
            if (args[i].indexOf("-tm") != -1)
            { // esli vremya
              i++;
              String s = args[i].trim();
              Double T = new Double(s);
              rna.Teta = T.doubleValue();
              continue;
            }
            if (args[i].indexOf("--time") != -1)
            { // esli vremya
              i++;
              String s = args[i].trim();
              Double T = new Double(s);
              rna.Teta = T.doubleValue();
              continue;
            }
            if (args[i].indexOf("-gr") != -1)
            { // esli konstanta rosta
              i++;
              String s = args[i].trim();
              Double K = new Double(s);
              rna.k_gr = K.doubleValue();
              continue;
            }
            if (args[i].indexOf("--growth") != -1)
            { // esli konstanta rosta
              i++;
              String s = args[i].trim();
              Double K = new Double(s);
              rna.k_gr = K.doubleValue();
              continue;
            }
            if (args[i].indexOf("-ncl") != -1)
            { // esli konstanta obrazovaniya pary
              i++;
              String s = args[i].trim();
              Double K = new Double(s);
              rna.k0 = K.doubleValue();
              continue;
            }
            if (args[i].indexOf("--nucleation") != -1)
            { // esli konstanta obrazovaniya pary
              i++;
              String s = args[i].trim();
              Double K = new Double(s);
              rna.k0 = K.doubleValue();
              continue;
            }
            if (args[i].indexOf("-m") != -1)
            { // esli kolichestvo prohodov
              i++;
              String s = args[i].trim();
              Integer I = new Integer(s);
              rna.M = I.intValue();
              continue;
            }
            if (args[i].indexOf("--m") != -1)
            { // esli kolichestvo prohodov
              i++;
              String s = args[i].trim();
              Integer I = new Integer(s);
              rna.M = I.intValue();
              continue;
            }
            if (args[i].indexOf("-out") != -1)
            { // esli imya vyhodnogo fila
              i++;
              String s = args[i].trim();
              rna.outputfile = s;
              continue;
            }
            if (args[i].indexOf("--outfile") != -1)
            { // esli imya vyhodnogo fila
              i++;
              String s = args[i].trim();
              rna.outputfile = s;
              continue;
            }
            if (args[i].indexOf("-hl") != -1)
            { // minimalnaya dlina spirali
              i++;
              String s = args[i].trim();
              Integer I = new Integer(s);
              rna.minhelix = I.intValue();
              continue;
            }
            if (args[i].indexOf("--helixLength") != -1)
            { // minimalnaya dlina spirali
              i++;
              String s = args[i].trim();
              Integer I = new Integer(s);
              rna.minhelix = I.intValue();
              continue;
            }
            if (args[i].indexOf("-gu") != -1)
            { // uchityvaem koncevye GU-pary
              i++;
              String s = args[i].trim();
              Integer B = new Integer(s);
              if (B.intValue() == 1) rna.if_GU = true;
              if (B.intValue() == 0) rna.if_GU = false;
              continue;
            }
            if (args[i].indexOf("--gu") != -1)
            { // uchityvaem koncevye GU-pary
              i++;
              String s = args[i].trim();
              Integer B = new Integer(s);
              if (B.intValue() == 1) rna.if_GU = true;
              if (B.intValue() == 0) rna.if_GU = false;
              continue;
            }
            if (args[i].indexOf("-ht") != -1)
            { // esli konstanta obrazovaniya pary
              i++;
              String s = args[i].trim();
              Double K = new Double(s);
              rna.helix_threshold = K.doubleValue();
              continue;
            }
            if (args[i].indexOf("--helixThreshold") != -1)
            { // esli konstanta obrazovaniya pary
              i++;
              String s = args[i].trim();
              Double K = new Double(s);
              rna.helix_threshold = K.doubleValue();
              continue;
            }

          }

          if (rna.elongate() == -1)
          {
            System.err.println("works too long");
          }
        }
        catch (OutOfMemoryError ex)
        {
          WriteFile writelog = new WriteFile(ID + ".log", false);
          writelog.Writeln("Out of memory");
        }
  /*      catch (Exception e)
        {
          WriteFile writelog = new WriteFile(ID + ".log", false);
          writelog.Writeln("Uncaught Exception");
        }//*/
      }

    //=================================================
      // VYDAET STRUKTURU LUCHSHEGO FOLDINGA
      foldings get_best_structure()
      {
        foldings f = null;
        group st1 = null;
        if (!st_time.isEmpty())
        {
          Set key_set = st_time.keySet(); // vytaskivaem set klyuchei
          Iterator iter = key_set.iterator();
          Float key = (Float) iter.next();
          st1 = (group) st_time.get(key);
          // esli prosto posledovatel'nost', to berem sleduyushii element
          if ((st1.elements[st1.best].folding_energy == 0)&&(iter.hasNext()))
          {
            key = (Float) iter.next();
            st1 = (group) st_time.get(key);
          }
        }
        f = st1.elements[st1.best];
        return f;
      }

      //=================================================
      // VYDAET STRUKTURU ZADANNOGO FOLDINGA
      group get_structure(int k)
      {
        if (k > st_time.size()) return null;
        group st1 = null;
        if (!st_time.isEmpty())
        {
          Set key_set = st_time.keySet(); // vytaskivaem set klyuchei
          Iterator iter = key_set.iterator();
          Float key = null;
          for (int i = 0; i <= k; i++)
            key = (Float) iter.next();
          st1 = (group) st_time.get(key);
        }
        return st1;
      }

      // vydaet element states, sootvetstvuyuschii dannomu exists
      foldings contain_structure(boolean[] exists)
      {
        group st1 = null;
        foldings f = null;
        if (groups_list1 != null)
        {
          if (groups_list1.hash.containsKey(exists))
          { // esli est', vytaskivaem nomer gruppy
            Integer group_index = (Integer) groups_list1.hash.get(exists);
            // i samu gruppu po etomu nomeru
            st1 = (group) groups_list1.get(group_index.intValue());
          }
        }
        if (st1 != null && st1.completed)
        f = st1.elements[st1.best];

        return f;
      }


      //======================================================================
      // vydaet element states, sootvetstvuyuschii pravilnoi strukture
      group find_by_pattern(String pattern)
      {
        Set key_set = st_time.keySet(); // vytaskivaem set klyuchei
        group st1 = null;
        if (!st_time.isEmpty())
        {
          for (Iterator iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
          {
            boolean[] key = (boolean[]) iter.next();
            st1 = (group) st_time.get(key);
            String str = st1.elements[st1.best].bracket_dots(sequence.length());
            if (str.compareTo(pattern) == 0)
            {
              return st1;
            }
          }
        }
        return null;
      }
      //==========================================================================
      // vydaet element states, kotoryi klevernyi list
       group clover()
       {
         Set key_set = st_energy.keySet(); // vytaskivaem set klyuchei
         group st1 = null;
         if (!st_time.isEmpty())
         {
           for (Iterator iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
           {
             Float key = (Float) iter.next();
             st1 = (group) st_energy.get(key);
             if (st1.elements[st1.best].clover_check())
             {
               return st1;
             }
           }
         }
         return null;
       }


    }
