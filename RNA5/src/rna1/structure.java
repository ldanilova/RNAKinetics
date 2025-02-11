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


// klass Comparator, po kotoromy budut sravnivat'sya sostoyaniya structury v TreeMap groups_list
class Comparator_exists implements Comparator
{
  public int compare (Object  a, Object b)
  {
    boolean[] a1 = (boolean[]) a;
    boolean[] b1 = (boolean[]) b;
//    for (int i = 0; i< a1.length; i++)
    for (int i = a1.length-1; i>=0; i--)
    {
      if (a1[i] == b1[i]) continue;
      if (b1[i])
      {// massiv a1 men'she
        return -1;
      }
      // massiv a1 bol'she
        return 1;
    }
    return 0;
  }
}
// Komporator dlya spiska petel' (loops_list)
// sravnivaet 2 vectora s nomerami spiralei. Vectory uporyadocheny gp vozrastaniyu
class Comparator_vector implements Comparator
{
  public int compare (Object  a, Object b)
  {
    Vector a1 = (Vector) a;
    Vector b1 = (Vector) b;
    if (a1.size()< b1.size()) return -1;
    if (a1.size()> b1.size()) return 1;
    for (int i = 0; i < a1.size(); i++)
    {
    	Integer ai = (Integer) a1.get(i);
    	Integer bi = (Integer) b1.get(i);
    	if (ai.compareTo(bi) == 0) continue;
        	return ai.compareTo(bi);
    }
    return 0;
  }
}

class Comparator_time implements Comparator
{ // po ubyvaniyu
  public int compare (Object  a, Object b)
  {
    double a1 = ((Float)a).floatValue();
    double b1 = ((Float)b).floatValue();
    if (a1 == b1) return 0;
    if (a1 > b1)
    {//  a1 bol'she
      return -1;
    }
    //  a1 men'she
    return 1;
    }
}

  class Comparator_energy implements Comparator
{ // po vozrastaniyu
  public int compare (Object  a, Object b)
  {
    double a1 = ((Float)a).floatValue();
    double b1 = ((Float)b).floatValue();
    if (a1 == b1) return 0;
    if (a1 < b1)
    {//  a1 men'she
      return -1;
    }
    //  a1 bol'she
    return 1;
    }
}


public class structure {
  final double k_gr = rna1.k_gr; // konstanta skorosti rosta cepi, t.e. poyavleniya ocherednogo nukleotida (c^-1)
  final double Teta = rna1.Teta; // process prekraschaetsya, kogda time stanet ravno Teta (c)

  // final int M = rna1.M; // kolichestvo strukturnyh perestanovok
  final int S = 500; // kolichestvo hranimyh sostoyanii, cherez kotorye my proshli (razmer TreeMap)
  final int min_stack = rna1.min_stack; // min kolichestvo sparennyh nukleotidov v strukture so skol'zyaschei petlei

  String ID = rna1.ID; // identificatior zadachi, imya log-file

  // describe RNA strand, possible helices and their dissipation constants
  heliset heliset1;
  int heliset_size;
  boolean exists[];
  int helix_rs2[];
  int RNA_length;
  double k_break[];

  foldings folding;
  double[][] graph_3;
  int ExistIndex = -1; // index v massivah exists, k_break, k_curr, s kotorymi mozhno rabotat' na dannyi moment

  // t.e. te spirali, kotorye na dannyi moment mogut suschestvovat'
  int CurPos = 0; // pokazyvaet do kakoi pozicii vyrosla cep'
  int GlobalCount = 0; // schitaet vse perestanovki
  int count_new = 0; // schitaet kolichestvo new v step

  double time = 0;

  //==========================================

  public structure(heliset in_heliset, int length)
  {
    RNA_length = length;
    heliset1 = in_heliset;
    heliset_size = heliset1.heliset_size;
    exists = new boolean[heliset_size];
    graph_3 = new double[rna1.graph_seq][length]; // graph_seq krivyh

    helix_rs2 = new int[heliset_size]; // massiv nachal pravogo plecha spirali


    // zapolnyaem massiv konechnyh pozicii spirali
    for (int i = 0; i < heliset_size; i++)
    {
      helix sh = (helix) heliset1.get(i);
      helix_rs2[i] = sh.D;
    }
  }

  //============================================
  // rostit RNK i skladyvaet ee v structury do vremeni Teta
  void run(groups_list groups_list1)
  {
    time = 0;
    CurPos = 0;
    ExistIndex = -1;
    // nachinaem s togo chto net niodnoi spirali
    for (int i = 0; i < heliset_size; i++) exists[i] = false;
    folding = new foldings(exists, heliset1);
    Vector exist_helises;

    //======================================
    // osnovnoi cikl po vremeni
    while (time < Teta)
    { // process prekraschaetsya, kogda time stanet ravno Teta
      // nahodim tekuschii ExistIndex
      for (int i = ExistIndex + 1; i < heliset_size; i++)
      {
        if (helix_rs2[i] <= CurPos)
        {
          ExistIndex = i;
          helix helix1 = (helix) heliset1.get(i);
          if (helix1.grown_time == 0)
          {
            helix1.grown_time = time;
          }
        }
      }
      //============================
      // zapominaem vremya, kogda vyrosla vsya posl-st'
      if ( (CurPos == RNA_length - 1) && (rna1.grown_time == 0))
      {
        rna1.grown_time = time;
      }
      //============================
      // shag v novoe sostoyanie, kotoryi proizoshel za vremya LiveTime
      double LiveTime = step(time, groups_list1);

      // poluchaem nomera suschestvuyuschih spiralei
      exist_helises = folding.index_exists_true();
      // kladem vremya v grafik dly spiralei (graph 2)
      if (!exist_helises.isEmpty()) heliset1.add_to_graph(time, LiveTime, exist_helises);

      double time_start = time; //nachalo zhizni struktury
      double time_end = time + LiveTime;
      if (time + LiveTime > Teta) time_end = Teta; // konec zhizni
      //============================================================
      // zapolnyaem grafik dlya posledovatel'nostei
      // v pozicii suschestvuyuschih spiralei pribavlyaem vremya
      // grafik 3: po x - posledovatel'nost', po y - sostoyanie dannoi pozicii (otkryta ili net)
      // dlya graph_seq momentov vremeni, t.e. za vremya [0, Teta/graph_seq], [Teta/graph_seq, 2*Teta/graph_seq]

      if (!exist_helises.isEmpty())
      { // esli v structure spirali est'
        int n_start = (int) ( (rna1.graph_seq * time_start) / Teta); // element massiva graph_3, kuda popadaet nachalo
        int n_end = (int) ( (rna1.graph_seq * time_end) / Teta); // element massiva graph_3, kuda popadaet konec
        for (int i = n_start; i <= n_end; i++)
        {
          double next_time = 0;
          if (i + 1 < rna1.graph_seq)
          {
            next_time = ( (i + 1) * Teta) / rna1.graph_seq; // vremya, sootvetstvuyuschee elementu i+1 v massive graph
          }
          else
          {
            for (int j = 0; j < exist_helises.size(); j++)
            {
              Integer I = (Integer) exist_helises.get(j);
              helix helix1 = (helix) heliset1.get(I.intValue());
              for (int k = helix1.A; k <= helix1.B; k++)
              { // levoe plecho
                if (time_end > Teta)
                  graph_3[rna1.graph_seq - 1][k] = graph_3[rna1.graph_seq - 1][k] + (Teta - time_start);
                else //time_end > Teta
                  graph_3[rna1.graph_seq - 1][k] = graph_3[rna1.graph_seq - 1][k] + (time_end - time_start);
              }
              for (int k = helix1.C; k <= helix1.D; k++)
              { // pravoe plecho
                if (time_end > Teta)
                  graph_3[rna1.graph_seq - 1][k] = graph_3[rna1.graph_seq - 1][k] + (Teta - time_start);
                else //time_end > Teta
                  graph_3[rna1.graph_seq - 1][k] = graph_3[rna1.graph_seq - 1][k] + (time_end - time_start);
              }
            }
            break;
          }
          if (next_time > time_end) { // esli konec zhizni lezhit v predydushem elemente
            // sootvetstvuyuschii otrezok vremeni dobavlyaem k elementu i
            for (int j = 0; j < exist_helises.size(); j++)
            {
              Integer I = (Integer) exist_helises.get(j);
              helix helix1 = (helix) heliset1.get(I.intValue());
              for (int k = helix1.A; k <= helix1.B; k++)
              { // levoe plecho
                graph_3[i][k] += (time_end - time_start);
              }
              for (int k = helix1.C; k <= helix1.D; k++)
              { // pravoe plecho
                graph_3[i][k] += (time_end - time_start);
              }
            }
          }
          else { // esli konec zhizni lezhit dal'she
            // sootvetstvuyuschii otrezok vremeni dobavlyaem k elementu i
            for (int j = 0; j < exist_helises.size(); j++)
            {
              Integer I = (Integer) exist_helises.get(j);
              helix helix1 = (helix) heliset1.get(I.intValue());
              for (int k = helix1.A; k <= helix1.B; k++)
              { // levoe plecho
                graph_3[i][k] += (next_time - time_start);
              }
              for (int k = helix1.C; k <= helix1.D; k++)
              { // pravoe plecho
                graph_3[i][k] += (next_time - time_start);
              }
            }
            // i sdvigaem vremya nachala
            time_start = next_time;
          }
        }
      }
      // sdvigaem tecuschee vremya na konec zhizni shpil'ki
      time = time_end;
      // uvelichivaem obschii schetchik perestanovok
      GlobalCount++;
    }
  }

  //==================================
  // shag v novoe sostoyanie
  double step(double CurrTime, groups_list groups_list1)
  {
    double k = 0; // parametr dlya Monte-Carlo = summe konstant obrazovaniya i raspada spiralei
    group gr = null; // element v folding_groups

    if (CurPos == 72)
    {
      CurPos = 72;
    }//*/

    // proveryaem net li gruppy s takim sostoyaniem
     if (groups_list1.hash.containsKey(folding.exists))
     {
       // esli est', vytaskivaem nomer gruppy
       Integer group_index = (Integer) groups_list1.hash.get(folding.exists);
       // i samu gruppu po etomu nomeru
       gr = (group) groups_list1.get(group_index.intValue());

  /*     if (folding.owner_group == null)
       {
         System.err.println("group # = " + group_index.intValue());
         System.err.println("folding " + "\n" + folding.toString(RNA_length));
       }//*/

       if (gr == null)
       {
         // pechataem group_list
         WriteFile out = new WriteFile("group_list.txt",false);
         out.Writeln(groups_list1.toString());
         // pechataem hash
         WriteFile out_hash = new WriteFile("hash.txt",false);
         Set key_set =  groups_list1.hash.keySet();
         Iterator iter;
         for (iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
         {
           boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
           Integer group_index1 = (Integer) groups_list1.hash.get(key_exists);
           StringBuffer sb = new StringBuffer();
           for (int j = 0; j < key_exists.length; j++)
           {
             if (key_exists[j]) sb.append("t" + "\t");
             else sb.append("f" + "\t");
           }
           sb.append("\n" + "Group# = " + group_index1.intValue());

           out_hash.Writeln(sb.toString());
         }
       }
       // proveryaem sobrana eta gruppa ili net
       if (gr != null && !gr.completed)
       {
         // esli net, to udalyaem ego iz spiska
         groups_list1.remove_group(group_index.intValue());
         // sobiraem
         gr.complete(folding);
         // i kladem obratno
         groups_list1.add_new_group(group_index.intValue(),gr);
  /*       rna1.writestruct.Writeln("GROUP");
         rna1.writestruct.Writeln(gr.toString());//*/
       }


     // uvelichivaem schetchik gruppy
      if (gr != null)
      {
        gr.count++; // uvelichivaem u nego schetchik
        gr.step = GlobalCount; // zapominaem novyi shag
      }
    }
    //==================================
    // esli takogo elementa eshe ne bylo
    else
    {
 //     long t1 = System.currentTimeMillis();

      // sozdaem group ot etogo foldinga
      gr = new group(folding, true);

 //     long t2 = System.currentTimeMillis();

      // dobavlyaem ego v spisok
      groups_list1.add_new_group(gr);
      count_new++;
    }


//    long t1 = System.currentTimeMillis();
    // nahodim konstantu dlya Monte-Carlo i perehod ili naraschivaem cep'
    Vector v = gr.get_perehod(ExistIndex,CurPos);
//    long t2 = System.currentTimeMillis();
    // vytaschili konstantu
    k = ((Double)v.get(0)).doubleValue();
    if (k == 0) return Teta;
    // vytaschili perehod
    foldings perehod = (foldings) v.get(1);
    CurPos = ((Integer)v.get(2)).intValue();
    // esli perehod nepustoi, to my ego vybiraem
    // esli pustoi, znachit my prosto narostili cep'
    if (perehod != null)
    {
      folding = perehod;
/*      //======================================
    StringBuffer sb = new StringBuffer();
      // zapisyvaem perehody

      for (int i = 0; i < heliset_size; i++)
      {
        sb.append(i + "\t");
      }
      sb.append("\n");

      // pechataem exists
      for (int i = 0; i < heliset_size; i++)
      {
        if (st.folding.exists[i]) sb.append("t" + "\t");
        else sb.append("f" + "\t");
      }
      sb.append("\n");
      // pechataem konstanty perehoda
      for (int i = 0; i < st.konst.length; i++)
      {
        sb.append(st.konst[i] + "\t");
      }
      sb.append("\n");
      sb.append("k = " + k + "\n" +
                "time = " + time + "\n" +
                "curr time = " + CurrTime + "\n" +
                "global count = " + GlobalCount + "\n" +
                folding.toString(RNA_length));


//    sb.append(folding.bracket_dots(RNA_length - 1));


      rna1.writestruct.Writeln(sb.toString());//*/


    }

    // vychislyaem vremya dannogo perehoda
    ExpRandom ExpRandom1 = new ExpRandom();
    double time = ExpRandom1.draw(k);
   // if (time < 0) System.err.println(time);
 //   if (time > 100) System.err.println(time);
//    if (CurPos >= RNA_length - 1) // budem sobirat' statistiku, esli cep' vyrosla
//    if (CurrTime >= 3) // budem sobirat' statistiku posle 3 sek
//    {
      if ( (CurrTime + time) > Teta) gr.live_time += (Teta - CurrTime);
      else gr.live_time += time; // uvelichivaem vremya
      gr.add_to_graph(CurrTime, time); // zapolnyaem grafik dlya structury
//    }
    return time; // vremya perehoda*/
//  return Teta;
  }
}
