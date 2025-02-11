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

class perehod
{
  foldings structure; // prehod v kakuy strukturu
  double constant; // s kakoi konstantoi

  public perehod(foldings f, double c)
  {
    structure = f;
    structure.down_loop = null;
    structure.up_loop = null;
    constant = c;
  }
}


/*class group_element
{
  final double R = 1.987; //gazovaya postoyannaya 8.314 Dj/K*mol' (1.987 kal/K*mol')
  final double T = 310; // temperature v K
  final double RT = R * T;

  foldings base_folding; // folding, odin element gruppy
  foldings[] perehody; // vozmozhnye perehody iz base_folding naruzhu
  double[] outside_const; // konstanty etih perehodov
  double weight; // ves etogo foldinga v gruppe
  boolean checked; // prosmotren etot element ili net
  int group_index; // pomnit iz kakoi gruppy etot element

  public group_element(foldings f)
  {
    base_folding = f;
    perehody = new foldings[base_folding.heliset_size];
    outside_const = new double[base_folding.heliset_size];
    weight = Math.exp(-base_folding.folding_energy/RT);
    checked = false;
    group_index = -1;
  }

  public String toString()
  {
    StringBuffer sb = new StringBuffer();
    // tochecno-skobochnaya zapis' struktury
    sb.append("weight = " + weight + "\n" + base_folding.toString(rna1.sequence.length()));
    // spirali iz etoi struktury
    sb.append("\n" + base_folding.print_helises(rna1.print_only_helix_index));
    return sb.toString();
  }
}//*/


class group
{
  final double k0 = rna1.k0; // konstanta skorosti obrazovaniya ocherednoi pary (c^-1)
  final double R = 1.987; //gazovaya postoyannaya 8.314 Dj/K*mol' (1.987 kal/K*mol')
  final double T = 310; // temperature v K
  final double RT = R * T;
  final double k_gr = rna1.k_gr; // konstanta skorosti rosta cepi, t.e. poyavleniya ocherednogo nukleotida (c^-1)
  double threshold_dis = rna1.threshold_dis; // porog na konstantu razrusheniya, dlya obrazovaniya gruppy
  double threshold_form = rna1.threshold_form; // porog na konstantu obrazovaniya

  foldings[] elements; // okonchatel'nyi nabor foldingov, kotorye vhodyat v odnu gruppu,
  perehod[] perehody; // perehody iz gruppy s konstantami
  double[] weights; // vesa elementov v gruppe
  private TreeMap processing_perehody; // promezhutochnyi nabor perehodov

  // poluchennyi iz TreeMap processing_elements
  TreeMap processing_elements; // poka sobiraem gruppu, zapominaem elementy v TreeMap
  private int elements_length; // dlina massiva s elementami

  // foldinga u kotorogo v dannoi pozicii v exists stoit samyi pravyi true
  // esli v etoi pozicii net true, to budet stoyat' -1
  private int[] rigth_true_perehod; // massiv po dline massiva perehody, ukazyvaet samyi bol'shoi nomer v elements
  int best; // nomer luchshego predstavitelya gruppy v elements
  boolean completed = false; //gruppa sobrana polnost'yu ili net
  int heliset_size;
  heliset heliset1;
  // statistika
  int count = 0; // skol'ko raz pobyvali v etoi gruppe
  int step = 0; // nomer shaga, v kotorom eta struktura obrazovalas' poslednii raz
  float live_time = 0; // summarnoe vremya, skol'ko posledovatelnost' byla v etom sostoyanii
  double graph_1[]; // vremena poyavlenii structury
  int group_index = -1; // nomer v group_list

  public group(foldings f, boolean comp)
  {
    heliset_size = f.heliset_size;
    heliset1 = f.heliset1;
    Comparator_exists C_exists = new Comparator_exists(); // klyuch, po kotoromu budut uporyadocheny elementy hash
    processing_elements = new TreeMap(C_exists);

    if (!comp)
    { // sozdaem gruppu iz odnogo elementa, ne sobiray
      add_element_to_group(f);
    }
    else
    { // sobiraet gruppu, schitaet vse konstanty perehodov i sami perehody, vesa foldingov
      rigth_true_perehod = new int[heliset_size];
      for (int j = 0; j < heliset_size; j++)
      { // delaem vse elementy -1
        rigth_true_perehod[j] = -1;
      }
      complete(f);
//      f.heliset1 = null;
    }
  }

  void add_element_to_group(foldings f)
  {
    processing_elements.put(f.exists, f);
    f.owner_group = this;
  }

  // zapolnyaen gruppu
  void complete(foldings f)
  {
    collect_group(f);
    if (!processing_elements.isEmpty())
    {
      // prosmatrivaem elements, poka ne soberem vsyu gruppu, t.e. poka check_results ne budet true
      boolean check_result = check_elements();
      while (!check_result)
      {
        check_result = check_elements();
      }
      // govorim, chto gruppa sobrana
      completed = true;

      // delaem iz processing_elements massiv elements
      elements_length = processing_elements.size();
      elements = new foldings[elements_length];
      Set key_set = processing_elements.keySet(); // vytaskivaem set klyuchei
      Iterator iter;
      int i = 0;
      for (iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
      {
        boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
        // prosmatrivaem net li u nego v massive perehody, structur iz gruppy
        foldings el = (foldings) processing_elements.get(key_exists);
        if (!el.checked) System.err.println(
            "non-checked element in the completed group");
        if (el.owner_group == null) System.err.println(
            "owner_group is not spesified (1)");
        if (el.owner_group != this) System.err.println(
            "owner_group does not equal curr group");
        for (int j = 0; j < heliset_size; j++)
        { // esli struktura iz gruppy schitaetsya perehodom naruzh, to udalaem
          if (el.perehody[j] != null &&
              processing_elements.containsKey(el.perehody[j].exists))
          {
            boolean[] j_exists = el.perehody[j].exists;
            if (el.perehody[j].owner_group == null)
              System.err.println("owner_group is null");

              // esli takaya structura byla v spiske, to my ee udalyaem
            if (rna1.groups_list1.hash.containsKey(j_exists))
            {
              // esli est', vytaskivaem gruppy
              group gr1 = (group) el.perehody[j].owner_group;
              /* if (el.perehody[j].owner_group == this)
               System.err.println("owner_group of perehod is this group");//*/
              Integer group_index = (Integer) rna1.groups_list1.hash.get(
                  j_exists);
              if (gr1 != null && gr1 != this && gr1.completed)
              {
                System.err.println("Try to remove completed group (1) = " +
                                   group_index.intValue());
                System.err.println(el.toString());
              }
              // udalyaem etu gruppu iz spiska
              rna1.groups_list1.remove_group(group_index.intValue());
              el.perehody[j].owner_group = this;
//              System.err.println("Udalyaem pri prosmotre perehodov gruppu #" + group_index.intValue());
              //             System.err.println(gr1.toString());
            }
            el.perehody[j] = null;
            el.outside_const[j] = 0;
          }
        }
        // i prisvaevaem ocherednomu elementu massiva
        elements[i] = el;
        i++;
      }
      // i pereschityvaem vesa
      weights = new double[elements_length];
      double weight_summ = 0; // summa vseh vesov etoi gruppy
      // i nahodim nomer luchshego predstavitelya
      double best_energy = Double.MAX_VALUE;
      for (int j = 0; j < elements_length; j++)
      {
        weights[j] = Math.exp(elements[j].folding_energy / RT);
        weight_summ += weights[j];
        if (elements[j].folding_energy < best_energy)
        {
          best_energy = elements[j].folding_energy;
          best = j;
        }
      }
      // normiruem vesa u vseh elementov gruppy
      for (int j = 0; j < elements_length; j++)
      {
        weights[j] = weights[j] / weight_summ;
      }

    }
    else System.err.println("elements empty!");
    // zapolnayaem massiv perehodov
    // snachala vhe perehody skladyvaem v TreeMap
    // potom iz nego delaem massiv perehody
    Comparator_exists C_exists = new Comparator_exists(); // klyuch, po kotoromu budut uporyadocheny elementy hash
    processing_perehody = new TreeMap(C_exists);
    // zapolnyaem processing_perehod
    // i ischem max konstantu
    double k_max = 0;
    for (int i = 0; i < elements_length; i++)
    {
      for (int j = 0; j < heliset_size; j++)
      {
        if (elements[i].perehody[j] != null)
        {
          // sozdaem perehod - struktura i konstanta perehoda * na ves etogo perehoda v gruppe
          foldings el = elements[i].perehody[j];

          // if (!el.checked) System.err.println("non-checked element in the completed group" );
          if (el.owner_group == null)
          {
            System.err.println("owner_group is not specified (2)");
          }
          perehod p = new perehod(el, (elements[i].outside_const[j] * weights[i]));
          processing_perehody.put(el.exists, p);
          if (k_max < p.constant) k_max = p.constant;
        }
      }
      elements[i].perehody = null;
      elements[i].outside_const = null;
      elements[i].folding = null;
      elements[i].up_loop = null;
      elements[i].down_loop = null;
     }

    rigth_true_perehod = new int[heliset_size];
    for (int j = 0; j < heliset_size; j++)
    { // delaem vse elementy -1
      rigth_true_perehod[j] = -1;
    }//*/

    if (k_max == 0) System.err.println("k_max = 0");
    Vector v = new Vector();
    if (!processing_perehody.isEmpty())
    {
      // ubiraem perehody s malen'koi konstantoi
      Set key_set = processing_perehody.keySet(); // vytaskivaem set klyuchei
      Iterator iter;
      for (iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
      {
        boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
        perehod pr = (perehod) processing_perehody.get(key_exists);
        // esli konstanta ne ochen' malen'kay, to kladem v perehody
        if (pr.constant / k_max > rna1.group_perehod)

          //processing_perehody.remove(key_exists);
          v.add(pr);
      }

      if (!v.isEmpty())
      { // delaem massiv perehody
        int size = v.size();
        perehody = new perehod[size];
   /*     for (int i = 0; i < size; i++)
        {
          perehod p = (perehod) v.get(i);
          perehody[i] = p;
        }//*/
        perehody = (perehod[])v.toArray(perehody);
        v = null;

        // zapolmyaem massiv rigth_true_perehod

        // prosmatrivaem exists s konca u kazhdogo elementa i zapominaem pervyi true
        for (int i = 0; i < size; i++)
          for (int j = heliset_size - 1; j >= 0; j--)
          {
            if (!perehody[i].structure.exists[j])continue;
            else
            {
              rigth_true_perehod[j] = i;
              break;
            }
          }
      }
    }

    processing_perehody = null;
    processing_elements = null;
    /*     System.err.println("gc in complete");
         System.gc();//*/
    graph_1 = new double[rna1.graph_step];
   for (int i = 0; i < rna1.graph_step; i++)
   {
     graph_1[i] = 0;
   }

  }

  //=================================================
  // vydaet vektor, v 0-om elemente kotorogo budet konstanta k dlya ExpRandom
  // v 1-om perehod (ot indeksa massiva exists i tekushei pozicii v cepi)
  // vo 2-om tekuschaya dlina cepi CurrPos
  /*  Vector get_perehod(int exist_index, int curr_pos)
    {
      Vector result = new Vector();
      double k = 0; // vozvraschaemaya konstanta
//   double k1 = 0; // konstanta dlya ArrayRandom
      foldings perehod = null; // vozvraschaemyi perehod
      int i_perehoda = -1; // nomer perehoda
      ArrayRandom ArrayRandom1 = new ArrayRandom(); // opredelyaet nomer sleduyuschego perehoda
      double[] curr_consts = null; // konstanty dlya k
      int curr_element = -1; // govorit, do kakogo indeksa mozhno brat' elementy iz massiva elements

      // hahodim nomer structury, kotoraya mozhet slozhit'sya na tekuschii moment
      // esli exist_index = -1, znachit esche net niodnoi spirali
       // nado esche rasti
       if (exist_index == -1)
       {
         result.add(0, new Double(k_gr));
         result.add(1, null);
         curr_pos++;
         result.add(2, new Integer(curr_pos));
         return result;
       }
   if (rigth_true[exist_index] != -1) curr_element = rigth_true[exist_index];
       else // ischem blizhaishii element ne ravnyi -1
         for (int i = exist_index - 1; i >= 0; i--)
           if (rigth_true[i] != -1)
           {
             curr_element = rigth_true[i];
             break;
           }
           else continue;

       if (curr_element == -1) // znachit maximal'no dopustimaya struktura - eto prosto posledovatel'nost'
       { // i nado esche rasti
         result.add(0, new Double(k_gr));
         result.add(1, null);
         curr_pos++;
         result.add(2, new Integer(curr_pos));
         return result;
       }

     // zadaem razmer massiva konstant, kak (exist_index*curr_element)
      if (exist_index < heliset_size - 1) // esli cep' esche ne vyrosla i nuzhna konstanta rosta
      {
        curr_consts = new double[ (exist_index+1) * (curr_element + 1) + 1];
        // i poslednii element = konstante rosta
        curr_consts[((exist_index+1) * (curr_element + 1))] = k_gr;
      }
      else // esli cep' vyrosla
        curr_consts = new double[ (exist_index+1) * (curr_element + 1)];
      // zapolnyaem massiv konstant

      for (int j = 0; j<= curr_element; j++ )// perebiraem massiv elements
      {
        for (int i = 0; i <= exist_index; i++) // perebiraem massiv konstant
          curr_consts[(exist_index+1)*j + i] = elements[j].outside_const[i]*elements[j].weight;
        j++;
      }
      // schitaem kostantu dlya Monte-Carlo
      for (int i = 0; i < curr_consts.length; i++)
      {
        k = k + curr_consts[i];
      }
      // ischem nomer perehoda
      i_perehoda = ArrayRandom1.draw(curr_consts, k);
      if (i_perehoda == -1)
      {
        WriteFile writelog; // log file
        writelog = new WriteFile(rna1.ID + ".log", true);
        writelog.Writeln("Number of next step doesn't find");
        System.err.println("Number of next step doesn't find");
      }
      // esli nado narostit' cep'
      if (i_perehoda == ((exist_index+1) * (curr_element + 1))) curr_pos++;
      else
      {
        // j_curr_elements - nomer v elements
        int j_curr_elements = (int)  i_perehoda/(exist_index+1);
        // ostatok ot deleniya = indeksu v massive perehodov
        int index_perehoda = i_perehoda % (exist_index + 1);
        perehod = elements[j_curr_elements].perehody[index_perehoda]; // berem sootvetstuyuschii perehod

      }
   /*    if (perehod != null){
          //======================================
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
            if (perehod.exists[i]) sb.append("t" + "\t");
            else sb.append("f" + "\t");
          }
          sb.append("\n");
          // pechataem konstanty perehoda
          for (int i = 0; i < curr_consts.length; i++)
          {
            sb.append(curr_consts[i] + "\t");
          }
          sb.append("\n");
          sb.append("k = " + k + " perehod = " + i_perehoda + "\n" +
                    perehod.toString(rna1.sequence.length()));

          rna1.writestruct.Writeln(sb.toString());
        }//


        // zapisyvaem v result vozvraschaemye znacheniya
        result.add(0, new Double(k));
        result.add(1, perehod);
        result.add(2, new Integer(curr_pos));
        return result;
      }//*/



    //=================================================
    // vydaet vektor, v 0-om elemente kotorogo budet konstanta k dlya ExpRandom
    // v 1-om perehod (ot indeksa massiva exists i tekushei pozicii v cepi)
    // vo 2-om tekuschaya dlina cepi CurrPos
    Vector get_perehod(int exist_index, int curr_pos)
    {
      Vector result = new Vector();
      double k = 0; // vozvraschaemaya konstanta
      foldings next_perehod = null; // vozvraschaemyi perehod
      int i_perehoda = -1; // nomer perehoda
      ArrayRandom ArrayRandom1 = new ArrayRandom(); // opredelyaet nomer sleduyuschego perehoda
      double[] curr_consts = null; // konstanty dlya k
      int curr_element = -1; // govorit, do kakogo indeksa mozhno brat' elementy iz massiva elements

      // hahodim nomer structury, kotoraya mozhet slozhit'sya na tekuschii moment
      // esli exist_index = -1, znachit esche net niodnoi spirali
      // nado esche rasti
      if (exist_index == -1)
      {
        result.add(0, new Double(k_gr));
        result.add(1, null);
        curr_pos++;
        result.add(2, new Integer(curr_pos));
        return result;
      }
      if (rigth_true_perehod[exist_index] != -1) curr_element =
          rigth_true_perehod[exist_index];
      else // ischem blizhaishii element ne ravnyi -1
        for (int i = exist_index - 1; i >= 0; i--)
          if (rigth_true_perehod[i] != -1)
          {
            curr_element = rigth_true_perehod[i];
            break;
          }
          else continue;

      if (curr_element == -1) // znachit maximal'no dopustimaya struktura - eto prosto posledovatel'nost'
      { // i nado esche rasti
        result.add(0, new Double(k_gr));
        result.add(1, null);
        curr_pos++;
        result.add(2, new Integer(curr_pos));
        return result;
      }

      // zadaem razmer massiva konstant
      if (exist_index < heliset_size - 1) // esli cep' esche ne vyrosla i nuzhna konstanta rosta
      {
        curr_consts = new double[curr_element + 2];
        // i poslednii element = konstante rosta
        curr_consts[curr_element + 1] = k_gr;
      }
      else // esli cep' vyrosla
        curr_consts = new double[perehody.length];
        // zapolnyaem massiv konstant
      for (int i = 0; i <= curr_element; i++)
        curr_consts[i] = perehody[i].constant;

        // schitaem kostantu dlya Monte-Carlo
      for (int i = 0; i < curr_consts.length; i++)
      {
        k = k + curr_consts[i];
      }
      // ischem nomer perehoda
      i_perehoda = ArrayRandom1.draw(curr_consts, k);
      if (i_perehoda == -1)
      {
        WriteFile writelog; // log file
        writelog = new WriteFile(rna1.ID + ".log", true);
        writelog.Writeln("Number of next step doesn't find");
        System.err.println("Number of next step doesn't find");
      }
      // esli nado narostit' cep'
      if (i_perehoda == curr_element + 1) curr_pos++;
      else
      {
        // sozdaem folding
        next_perehod = perehody[i_perehoda].structure; // berem sootvetstuyuschii perehod
      }
      /*    if (perehod != null){
            //======================================
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
              if (perehod.exists[i]) sb.append("t" + "\t");
              else sb.append("f" + "\t");
            }
            sb.append("\n");
            // pechataem konstanty perehoda
            for (int i = 0; i < curr_consts.length; i++)
            {
              sb.append(curr_consts[i] + "\t");
            }
            sb.append("\n");
            sb.append("k = " + k + " perehod = " + i_perehoda + "\n" +
                      perehod.toString(rna1.sequence.length()));

            rna1.writestruct.Writeln(sb.toString());
          }//*/


       // zapisyvaem v result vozvraschaemye znacheniya
       result.add(0, new Double(k));
      result.add(1, next_perehod);
      result.add(2, new Integer(curr_pos));
      return result;
    } //*/

    //===========================================
    // proveryaet u vseh li elementov zapolneny massivy perehodov i konstant
  boolean check_elements()
  {
    Set key_set = processing_elements.keySet(); // vytaskivaem set klyuchei
    for (Iterator iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
    {
      boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
      foldings gr_el = (foldings) processing_elements.get(key_exists); // izvlekaem sootvetstvuyuschii element
      // esli on ne zapolnen
      if (!gr_el.checked)
      { // to vyzyvaem collect ot etogo elementa
        collect_group(gr_el);
        return false;
      }
    }
    /*    System.err.println("gc in check_elements");
        System.gc();//*/
     return true;
  }

  // ==========================================================
  // nahodit perehod
  foldings find_perehod(foldings f)
  {
    foldings gr_el = null; // vozvraschaemyi perehod
    // esli gruppa s takim elementom uzhe est'
    // to pisvaevaem ee perehodu
    // vytaskivaem nomer gruppy
    Integer group_index = (Integer) rna1.groups_list1.hash.get(f.exists);
    group gr = (group) rna1.groups_list1.get(group_index.intValue());
    if (gr.completed)
    {
      // ischem v elements folding ravnyi folding_next
      Comparator_exists C_exists = new Comparator_exists(); //
      for (int i = 0; i < gr.elements.length; i++)
        if (C_exists.compare(f.exists, gr.elements[i].exists) == 0)

          // i prisvaivaem ego perehodu
          gr_el = gr.elements[i];
    }
    else
    {
      gr_el = (foldings) gr.processing_elements.get(
          f.exists);
    }
    if (gr_el == null) System.err.println("null element");

    return gr_el;
  }

//=============================================================
// sobiraet gruppu vokrug odnogo foldinga
  void collect_group(foldings current_folding)
  {
    double konst_dis, konst_form;
    // esli net niodnogo elementa ili est' no ne current_f
    if ( (!processing_elements.isEmpty() &&
          (!processing_elements.containsKey(current_folding.exists)))
        || processing_elements.isEmpty())
    {
      current_folding.checked = true; // govorim, chto my ego prosmotreli
      // i kladem current_folding kak pervyi element
      this.add_element_to_group(current_folding);
    }
    else // ili takoi element est',
    if (processing_elements.containsKey(current_folding.exists))
    {
      // esli zapolnen, to nichego ne delaem
      if (current_folding.checked)
      {
        return;
      }

    }
    // elsi ne zapolnen, to zapolnayem
    // perebiraem vse spirali i skladyvaem foldingi
    for (int t = 0; t < heliset_size; t++)
    {
      foldings folding_next = new foldings(current_folding);
      if (current_folding.perehody == null)
      {
        current_folding.perehody = new foldings[heliset_size];;
        current_folding.outside_const = new double[heliset_size];

      }
      helix helix2 = (helix) heliset1.get(t); // tekuschaya spiral
      if (current_folding.exists[t] == false)
      { // esli eta spiral mozhet obrazovat'sya
        //================================================
        // polnaya sovmestimost' spiralei
        boolean can_zip = true; // govorit mozhet li t-aya spiral' obrazovat'sya, t.e. ona sovmestima so vsemi ostal'nymi
        int uncompatible = 0; // schitaet so skol'kimi ne sovmestima
        int sliding_ends = 0; // 0 - 4
        // mi sobiraem vsie _sushestvuiushie_
        // spirali, s kotorimi nasha kandidatnaia nahoditsia v "sliding-end" otnosheniah
        for (int j = 0; j < heliset_size; j++)
        {
          int Compat = heliset1.Compatible[t][j];
          if (Compat == 0)continue;
          if ( (current_folding.exists[j]) && (Compat == -1))
          {
            can_zip = false; // esli ne sovmestima hotya by s odnoi iz uzhe suschestvuyuschih
            uncompatible++;
          }
          if ( (current_folding.exists[j]) && (Compat > 0))
          {
            sliding_ends++;
          }
        }
        if (can_zip)
        { // esli t-aya spiral mozhet obrazovat'sya
          int n = helix2.Left.length(); // dlina tekuschei spirali

          // esli mogut obrazovat'sya sk. petli, to vyzyvaem funkciyu podscheta dliny
          if (sliding_ends > 0) n = heliset1.get_length_without_slid(t,
              current_folding.exists);
          // esli dlina kuska spirali bez sk. petlei >1, to takuyu spiral' obrazovyvaem
          if (n > 1)
          {
            folding_next.add_helix(t); // dobavlyaem spiral t
            // schitaem konstantu skorosti obrazovaniya (zavisit tol'ko ot energii petel' i dliny spirali!!!!)
            konst_form = k0 * (n - 1) *
                Math.exp( (current_folding.loops_energy -
                           folding_next.loops_energy) / (RT));
            // schitaem obratnuyu konstantu razrushtniya
            double[] get_energy_length = heliset1.
                get_length_energi_without_slid(t, folding_next.exists);
            double split_energy = get_energy_length[1];
            n = (int) get_energy_length[0];
            konst_dis = 0;
            if (n < 1 || split_energy > 0)
            {
              System.err.println("n = " + n + "\n" + "energy = " +
                                 split_energy
                                 /*+ "\n" + "t = " + t + "\n" + "global count = " + GlobalCount*/
                                 );
            }
            else
            {
              konst_dis = (n - 1) * k0 * Math.exp(split_energy / (RT));

              // sravnivaem s porogom obe konstanty, esli oni obe bol'she poroga
              // to eto perehod vnutri gruppy
              if ( (konst_dis > threshold_dis) && (konst_form > threshold_form))
              {
                // esli ego tam esche net
                if (!processing_elements.containsKey(folding_next.exists))
                {
                  // esli takaya structura byla v spiske, to my ee udalyaem
                  if (rna1.groups_list1.hash.containsKey(folding_next.
                      exists))
                  {
                    // esli est', vytaskivaem nomer gruppy
                    Integer group_index = (Integer) rna1.groups_list1.hash.
                        get(folding_next.exists);
                    group gr = (group) rna1.groups_list1.get(group_index.
                        intValue());
                    if (gr != null && gr.completed)
                    {
                      System.err.println("Try to remove completed group (2)" +
                                         group_index.intValue());
                      System.err.println(gr.toString());
                    }
                    // to kladem v elements tot folding, kotoryi byl v toi nezapolnennoi gruppe
                    // t.e. prosto perepisyvaem owner_group elementu gr
                    foldings gr_el = (foldings) gr.processing_elements.get(
                        folding_next.exists);
                    if (gr_el == null) System.err.println("null element");
                    else
                      add_element_to_group(gr_el);
                      // udalyaem etu gruppu iz spiska
                    rna1.groups_list1.remove_group(group_index.intValue());
                    folding_next = null;
                  }
                  else add_element_to_group(folding_next);
                  continue;
                }
              }
              else // esli hotya by odna men'she, to eto perehod vne gruppy
              {
                // elsi etu strukturu ecshe ran'she ne polozhili vnutr' tekuschei gruppy
                if (!processing_elements.containsKey(folding_next.exists))
                {
                  // zapominaem konstanru perehoda
                  current_folding.outside_const[t] = konst_form;
                  // i dobavlyaem folding_next v spisok grupp, kak nezapolnennyi element
                  // esli ego esche tam net
                  if (!rna1.groups_list1.hash.containsKey(folding_next.
                      exists))
                  { // my zapominaem etot perehod i ego konstantu
                    current_folding.perehody[t] = folding_next;
                    group gr1 = new group(folding_next, false);
                    rna1.groups_list1.add_new_group(gr1);
                  }
                  else
                  { // ischem perehod
                    current_folding.perehody[t] = this.find_perehod(
                        folding_next);
                  }
                }
                continue;
              }
            }
          }
          folding_next = null;
        }
        else
        {
          //  огда процесса два последовательно, складываютс€ не константы, а их  обратные величины (характерные времена)
          { // t-aya spiral' hot' s kem-to nesovmestima (can_zip = false)
            if (uncompatible > 0 /*&& sliding_ends <= 1*/)
            // esli nesovmestima s odnoi i net sk. kontsov, to ischem s kakoi
            {

              // dlya kostanty obrazovaniya (konst_form)
              double energy_break = 0;
              int n = 0; // kolichestvo par dlyz iniciacii

              // esli t-ya spiral' obrazuet s kem-to sk. petlyu
              if (sliding_ends > 0)
              {
                n = heliset1.get_length_without_slid(t,
                    current_folding.exists);
                // esli sk. petli s'edayut spiral'
                if (n <= 0)continue;
                n = 0;
              }

              // dlya konstanty razrusheniya (konst_dis)
              double energy_break_back = 0;
              int n_back = 0; // kolichestvo par dlya iniciacii obratnogo perehoda
              boolean find_const = false; // esli udalos' opredelit' tip nesovmestimosti
              // t.e. odin iz sluchaev podoshel
              double transit_energy = 0;
              int j1 = -1;
              // sozdaem massivy dlya perehodnogo sostoyaniya
              boolean[] t_exists = new boolean[heliset_size];
              pair[] pairsAD = new pair[heliset_size];
              pair[] pairsBC = new pair[heliset_size];
              // sozdaem massivy dlya obratnogo perehodnogo sostoyaniya
              boolean[] t_exists_back = new boolean[heliset_size];
              pair[] pairsAD_back = new pair[heliset_size];
              pair[] pairsBC_back = new pair[heliset_size];

              // energiya ostavshegosya kusochka spirali, esli vozmozhna iniciaciya pri ne polnom razrushenii
              double energy_of_part = 0;
              // zapolnyaem massivy dlya perehodnogo sostoyaniya
              for (int k1 = 0; k1 < heliset_size; k1++)
              {
                t_exists[k1] = current_folding.exists[k1];
                t_exists_back[k1] = current_folding.exists[k1];

                // zapolnyaem massivy pairsAD i pairsBC
                if (t_exists[k1])
                {
                  helix h = (helix) heliset1.get(k1);
                  pair p1 = new pair(h.A, h.D);
                  pair p2 = new pair(h.B, h.C);
                  pairsAD[k1] = p1;
                  pairsBC[k1] = p2;
                  pairsAD_back[k1] = p1;
                  pairsBC_back[k1] = p2;
                }
              }

              for (int j = 0; j < heliset_size; j++)
                if ( (current_folding.exists[j]) &&
                    (heliset1.Compatible[t][j] == -1))
                { // nashli spiral', kotoraya suschestvuet i nesovmestima s t-oi
                  // proveraym kak nesovmestima
                  helix helix1 = (helix) heliset1.get(j); // uzhe est' v strukture
 //                 int n1 = helix1.Left.length(); // dlina suschestvuyuschei spirali
 //                 int n2 = helix2.Left.length(); // dlina spirali, kotoraya hochet obrazovat'sya
                  int A1 = helix1.A; // nachalo levogo plecha
                  int A2 = helix2.A;
                  int B1 = helix1.B; // konec levogo plecha
                  int B2 = helix2.B;
                  int C1 = helix1.C; // nachalo pravogo plecha
                  int C2 = helix2.C;
                  int D1 = helix1.D; // konec pravogo plecha
                  int D2 = helix2.D;
                  String left;
                  String right;

                  //==================================
                  // SLUCHAI 5
                  // neobhodimo polnoe razrushenie suschestvuyuschei j-oi (sluchai 5)
                  if ( ( (B2 <= A1) && (D2 <= C1) && (C2 >= B1)) ||
                      ( (A2 >= B1) && (C2 >= D1) && (B2 <= C1)))
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    // esli est' sk. petli, to razrushat' mozhet byt' nado men'she
                    /*  if (sliding_ends > 0)
                      {
                        double[] get_energy_length = heliset1.
                            get_length_energi_without_slid(t,
                            current_folding.exists);
                        double split_energy = get_energy_length[1];
                         n = (int) get_energy_length[0];
                         energy_break += split_energy;
                      } else//*/
                    energy_break += helix1.energy;
                    // j-ya razrushina
                    t_exists[j] = false;
                    pairsAD[j] = null;
                    pairsBC[j] = null;
                    // ot t-oi odna para
                    // esli my esche ee ne prisvoili
                    if (!t_exists[t])
                    {
                      t_exists[t] = true;
                      pair p = new pair(A2, D2);
                      pairsAD[t] = p;
                      pairsBC[t] = p;
                      p = null;

                    }
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      energy_break_back = helix2.energy;
                      n_back = 1;
                      // dlya obratnogo perehodnogo sostoyaniya
                      // ot j-oi odna para
                      t_exists_back[j] = true;
                      pair p = new pair(A1, D1);
                      pairsAD_back[j] = p;
                      pairsBC_back[j] = p;
                      // t-aya razrushaetsya
                      t_exists_back[t] = false;
                      pairsAD_back[t] = null;
                      pairsBC_back[t] = null;
                      p = null;
                    }

                    find_const = true;
                    continue;
                  }

                  //============================================
                  //SLUCHAI 3
                  // vozmozhna svobodnaya iniciaciya i novaya spiral' vklyuchaet v sebya suschestvuyuschuyu
                  // novaya spiral' vytesnit staruyu (sluchai 3)
                  if ( (D2 >= D1) && (C2 < C1) && // pravoe plecho sp.1 vklucheno v pravoe plecho sp.2
                      (A2 - A1 > D1 - D2)) // naklon vpravo
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    // para, po kotoroi obrazuetsya perehodnoe sostoyanie
                    pair p = null;
                    pair p1_back = null;
                    pair p2_back = null;
                    pair p3_back = null;
                    //chislo par, godnih k initsiatsii spirali sverhu ot sush. 1-i (j-oi) spirali

                    n += Math.min( (B2 - B1), (C1 - C2));
                    p = new pair(B2, C2);
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      int n_break = D2 - D1 + 1; // chislo par, k-rye nado razrushit'
                      left = helix2.Left.substring(0, n_break);
                      right = helix2.Right.substring(0, n_break);
                      // energiya razrushennogo kusochka
                      energy_break_back = helix2.energy_calculation(left, right);

                      t_exists_back[j] = true;
                      p1_back = new pair(A1, D1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      p2_back = new pair(A2 + n_break, D2 - n_break);
                      p3_back = new pair(B2, C2);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;

                      p1_back = null;
                      p2_back = null;
                      p3_back = null;
                    }
                    if (n <= 0) System.err.println("n = " + n + "\t1");

                    // ot t-oi odna para
                    t_exists[t] = true;
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;

                    find_const = true;
                    continue;
                  }
                  if ( (D2 > D1) && (C2 <= C1) && // pravoe plecho sp.1 vklucheno v pravoe plecho sp.2
                      (A2 - A1 < D1 - D2)) // naklon vlevo
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    // para, po kotoroi obrazuetsya perehodnoe sostoyanie
                    pair p = null;
                    pair p1_back = null;
                    pair p2_back = null;
                    pair p3_back = null;

                    n += Math.min( (A1 - A2), (D2 - D1));
                    p = new pair(A2, D2);
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      // dlya obratnogo perehodnogo sostoyaniya
                      n_back = 1;
                      int n_break = C1 - C2 + 1; // chislo par, k-rye nado razrushit'
                      left = helix2.Left.substring(helix2.Left.length() -
                          n_break);
                      right = helix2.Right.substring(helix2.Left.length() -
                          n_break);
                      energy_break_back = helix2.energy_calculation(left,
                          right);

                      // ot j-oi odna para
                      t_exists_back[j] = true;
                      p1_back = new pair(B1, C1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;

                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      p2_back = new pair(A2, D2);
                      p3_back = new pair(B2 - n_break, C2 + n_break);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;

                    }
                    if (n <= 0) System.err.println("n = " + n + "\t2");

                    // ot t-oi odna para
                    t_exists[t] = true;
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;

                    find_const = true;
                    continue;
                  }

                  if ( (A2 < A1) && (B2 >= B1) && // levoe plecho sp.1 vklucheno v levoe plecho sp.2
                      (A2 - A1 > D1 - D2)) // naklon vpravo
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    pair p = null;
                    pair p1_back = null;
                    pair p2_back = null;
                    pair p3_back = null;
                    //chislo par, godnih k initsiatsii spirali snizu ot sush. 1-i (j-oi) spirali

                    n += Math.min( (A1 - A2), (D2 - D1));
                    p = new pair(A2, D2);
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      int n_break = B2 - B1 + 1; // chislo par, k-rye nado razrushit'
                      left = helix2.Left.substring(helix2.Left.length() -
                          n_break);
                      right = helix2.Right.substring(helix2.Left.length() -
                          n_break);
                      energy_break_back = helix2.energy_calculation(left, right);

                      t_exists_back[j] = true;
                      p1_back = new pair(B1, C1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      p2_back = new pair(A2, D2);
                      p3_back = new pair(B2 - n_break, C2 + n_break);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;
                    }

                    if (n <= 0) System.err.println("n = " + n + "\t3");

                    t_exists[t] = true;
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;

                    find_const = true;
                    continue;
                  } //

                  if ( (A2 <= A1) && (B2 > B1) && // levoe plecho sp.1 vklucheno v levoe plecho sp.2
                      (A2 - A1 < D1 - D2)) // naklon vlevo
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    pair p = null;
                    pair p1_back = null;
                    pair p2_back = null;
                    pair p3_back = null;

                    //chislo par, godnih k initsiatsii spirali sverhu ot sush. 1-i (j-oi) spirali

                    n += Math.min( (B2 - B1), (C1 - C2));
                    p = new pair(B2, C2);
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      int n_break = A1 - A2 + 1; // chislo par, k-rye nado razrushit'
                      left = helix2.Left.substring(0, n_break);
                      right = helix2.Right.substring(0, n_break);
                      // energiya razrushennogo kusochka
                      energy_break_back = helix2.energy_calculation(left, right);

                      t_exists_back[j] = true;
                      p1_back = new pair(A1, D1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      p2_back = new pair(A2 + n_break, D2 - n_break);
                      p3_back = new pair(B2, C2);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;
                    }

                    if (n <= 0) System.err.println("n = " + n + "\t4");

                    t_exists[t] = true;
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;

                    find_const = true;
                    continue;
                  } //

                  if ( (A2 <= C1) && (B2 > D1))
                  // pravoe plecho sp.1 vklucheno v levoe plecho sp.2
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    // para, po kotoroi obrazuetsya perehodnoe sostoyanie
                    pair p = null;
                    pair p1_back = null;
                    pair p2_back = null;
                    pair p3_back = null;
                    //chislo par, godnih k initsiatsii spirali snizu ot sush. 1-i (j-oi) spirali
                    n += (B2 - D1);
                    p = new pair(B2, C2);

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      //chislo par, kotorye nado razrushit
                      n_back = 1;
                      int n_break = C1 - A2 + 1; // chislo par, k-rye nado razrushit'
                      left = helix2.Left.substring(0, n_break);
                      right = helix2.Right.substring(0, n_break);
                      // energiya razrushennogo kusochka
                      energy_break_back = helix2.energy_calculation(left, right);

                      t_exists_back[j] = true;
                      p1_back = new pair(B1, C1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      p2_back = new pair(A2 + n_break, D2 - n_break);
                      p3_back = new pair(B2, C2);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;

                    }
                    if (n <= 0) System.err.println("n = " + n + "\t5");

                    t_exists[t] = true;
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;

                    find_const = true;
                    continue;
                  }

                  if ( (D2 >= B1) && (C2 < A1))
                  // levoe plecho sp.1 vklucheno v pravoe plecho sp.2
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    // para, po kotoroi obrazuetsya perehodnoe sostoyanie
                    pair p = null;
                    pair p1_back = null;
                    pair p2_back = null;
                    pair p3_back = null;

                    //chislo par, godnih k initsiatsii spirali snizu ot sush. 1-i (j-oi) spirali
                    n += (A1 - C2);
                    p = new pair(B2, C2);

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      int n_break = D2 - B1 + 1; // chislo par, k-rye nado razrushit'
                      left = helix2.Left.substring(0, n_break);
                      right = helix2.Right.substring(0, n_break);
                      // energiya razrushennogo kusochka
                      energy_break_back = helix2.energy_calculation(left, right);

                      t_exists_back[j] = true;
                      p1_back = new pair(B1, C1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      p2_back = new pair(A2 + n_break, D2 - n_break);
                      p3_back = new pair(B2, C2);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;

                    }

                    if (n <= 0) System.err.println("n = " + n + "\t6");

                    t_exists[t] = true;
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;

                    find_const = true;
                    continue;
                  }
                  //KONETS SLUCHAIA 3

                  //========================================
                  //OBRATNII SLUCHAI 3
                  // ne vozmozhna svobodnaya iniciaciya ni s odnogo kontsa
                  // i novaya spiral' vklyuchena v suschestvuyuschuyu
                  // novaya spiral' trebuet pokusa obeih kontsov;
                  // kinetiku takoi zarazi mi s hodu pridumat ne mojem, a potonu
                  // dalneishee est rezultat togo, cto mi verim v lokalnii balans otnositelno
                  // raspredelenia Boltzmana

                  if ( (D1 >= D2) && (C1 <= C2))
                  // pravoe plecho t-oi vklucheno v pravoe plecho j-oi
                  {
                    // sluchai 3A
                    if (A2 - A1 > D1 - D2) // naklon vpravo
                    {
                      folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                      int n_up = C2 - C1 + 1;
                      int n_down = A2 - A1 + 1;
                      if ( (n_down > 0 && n_up <= n_down) || (n_down <= 0))
                      {
                        //chislo par, kotorye nado razrushit
                        int n_break = n_up; // chislo par, k-rye nado razrushit'
                        left = helix1.Left.substring(helix1.Left.length() -
                            n_break);
                        right = helix1.Right.substring(helix1.Left.length() -
                            n_break);
                        energy_break += helix1.energy_calculation(left, right);

                        // energiya ostavshegosya kusochka
                        left = helix1.Left.substring(0,
                            helix1.Left.length() - n_break);
                        right = helix1.Right.substring(0,
                            helix1.Left.length() - n_break);
                        energy_of_part += helix1.energy_calculation(left, right);

                        t_exists[t] = true;
                        pair p = new pair(B2, C2);
                        pairsAD[t] = p;
                        pairsBC[t] = p;
                        // ot j-oi spirali est' kusochek
                        t_exists[j] = true;
                        pair p1 = new pair(A1, D1);
                        pair p2 = new pair(B1 - n_break, C1 + n_break);
                        pairsAD[j] = p1;
                        pairsBC[j] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }
                      if (n_down > 0 && n_up > n_down)
                      {
                        int n_break = n_down; // chislo par, k-rye nado razrushit' snizu
                        left = helix1.Left.substring(0, n_break);
                        right = helix1.Right.substring(0, n_break);
                        // energiya razrushennogo kusochka
                        energy_break += helix1.energy_calculation(left, right);

                        // energiya ostavshegosya kusochka
                        left = helix1.Left.substring(n_break);
                        right = helix1.Right.substring(n_break);
                        energy_of_part += helix1.energy_calculation(left, right);
                        t_exists[t] = true;
                        pair p = new pair(A2, D2);
                        pairsAD[t] = p;
                        pairsBC[t] = p;
                        // ot j-oi spirali est' kusochek
                        t_exists[j] = true;
                        pair p1 = new pair(A1 + n_break, D1 - n_break);
                        pair p2 = new pair(B1, C1);
                        pairsAD[j] = p1;
                        pairsBC[j] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }

                      // esli nesovmestima tol'ko s odnoi spiral'yu,
                      // to nam nuzhna budet konstanta obratnogo perehoda
                      if (uncompatible == 1)
                      {
                        if (D1 == D2 && B1 <= B2)
                        // znachit svobodnaya iniciaciya nevozmozhna
                        // nuzhno razrushit' odnu paru snizu
                        {
                          n_back = 1;
                          // energiya razrushennogo kusochka
                          energy_break_back = 0;
                          // ot j-oi odna para
                          pair p = new pair(A1, D1);
                          t_exists_back[j] = true;
                          pairsAD_back[j] = p;
                          pairsBC_back[j] = p;

                          t_exists_back[t] = true;
                          pair p1 = new pair(A2 + n_back, D2 - n_back);
                          pair p2 = new pair(B2, C2);
                          pairsAD_back[t] = p1;
                          pairsBC_back[t] = p2;
                          p = null;
                          p1 = null;
                          p2 = null;
                        }
                        else // kogda vozmozhna svobodnaya iniciaciya
                        {
                          n_down = 0;
                          pair p = null;
                          if (A2 > A1)
                          {
                            n_down = Math.min( (A2 - A1), (D1 - D2));
                            p = new pair(A1, D1);
                          }
                          n_up = 0;
                          if (B1 > B2)
                          {
                            n_up = Math.min( (B1 - B2), (C2 - C1));
                            p = new pair(B1, C1);
                          }
                          n_back = n_down + n_up;
                          if (n_back == 0) System.err.println("n = 0, obr. 3A");
                          t_exists_back[t] = true;
                          pair p1 = new pair(A2, D2);
                          pair p2 = new pair(B2, C2);
                          pairsAD_back[t] = p1;
                          pairsBC_back[t] = p2;

                          // ot j-oi odna para
                          t_exists_back[j] = true;
                          pairsAD_back[j] = p;
                          pairsBC_back[j] = p;
                          p = null;
                          p1 = null;
                          p2 = null;
                        }

                      }
                      find_const = true;
                      continue;
                    }
                    // sluchai 3B
                    if (A2 - A1 < D1 - D2) // naklon vlevo, razrushaem snizu
                    {
                      folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                      int n_down = D1 - D2 + 1;
                      int n_up = B1 - B2 + 1;
                      //chislo par, kotorye nado razrushit
                      if ( (n_up > 0 && n_down <= n_up) || (n_up <= 0))
                      {
                        int n_break = n_down; // chislo par, k-rye nado razrushit' snizu
                        left = helix1.Left.substring(0, n_break);
                        right = helix1.Right.substring(0, n_break);
                        // energiya razrushennogo kusochka
                        energy_break += helix1.energy_calculation(left, right);

                        // energiya ostavshegosya kusochka
                        left = helix1.Left.substring(n_break);
                        right = helix1.Right.substring(n_break);
                        energy_of_part += helix1.energy_calculation(left, right);
                        t_exists[t] = true;
                        pair p = new pair(A2, D2);
                        pairsAD[t] = p;
                        pairsBC[t] = p;
                        // ot j-oi spirali est' kusochek
                        t_exists[j] = true;
                        pair p1 = new pair(A1 + n_break, D1 - n_break);
                        pair p2 = new pair(B1, C1);
                        pairsAD[j] = p1;
                        pairsBC[j] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }
                      if (n_up > 0 && n_down > n_up)
                      { // razrushaem sverhu
                        int n_break = n_up;
                        left = helix1.Left.substring(helix1.Left.length() -
                            n_break);
                        right = helix1.Right.substring(helix1.Left.length() -
                            n_break);
                        energy_break += helix1.energy_calculation(left, right);
                        t_exists[t] = true;
                        pair p = new pair(B2, C2);
                        pairsAD[t] = p;
                        pairsBC[t] = p;
                        // ot j-oi spirali est' kusochek
                        t_exists[j] = true;
                        pair p1 = new pair(A1, D1);
                        pair p2 = new pair(B1 - n_break, C1 + n_break);
                        pairsAD[j] = p1;
                        pairsBC[j] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }

                      // esli nesovmestima tol'ko s odnoi spiral'yu,
                      // to nam nuzhna budet konstanta obratnogo perehoda
                      if (uncompatible == 1)
                      {
                        if ( (C1 == C2) && (A1 >= A2))
                        // znachit svobodnaya iniciaciya nevozmozhna
                        // nuzhno razrushit' odnu paru sverhu
                        {
                          n_back = 1;
                          // energiya razrushennogo kusochka
                          energy_break_back = 0;
                          pair p = new pair(B1, C1);
                          // ot j-oi odna para
                          t_exists_back[j] = true;
                          pairsAD_back[j] = p;
                          pairsBC_back[j] = p;

                          t_exists_back[t] = true;
                          pair p1 = new pair(A2, D2);
                          pair p2 = new pair(B2 - n_back, C2 + n_back);
                          pairsAD_back[t] = p1;
                          pairsBC_back[t] = p2;
                          p = null;
                          p1 = null;
                          p2 = null;
                        }
                        else // kogda vozmozhna svobodnaya iniciaciya
                        {
                          n_down = 0;
                          n_up = 0;
                          pair p = null;
                          if (A2 > A1)
                          {
                            n_down = Math.min( (A2 - A1), (D1 - D2));
                            p = new pair(A1, D1);
                          }
                          if (B1 > B2)
                          {
                            n_up = Math.min( (B1 - B2), (C2 - C1));
                            p = new pair(B1, C1);
                          }
                          t_exists_back[t] = true;
                          pair p1 = new pair(A2, D2);
                          pair p2 = new pair(B2, C2);
                          pairsAD_back[t] = p1;
                          pairsBC_back[t] = p2;
                          n_back = n_down + n_up;
                          if (n_back == 0) System.err.println("n = 0, obr. 3B");
                          // ot j-oi odna para
                          t_exists_back[j] = true;
                          pairsAD_back[j] = p;
                          pairsBC_back[j] = p;
                          p = null;
                          p1 = null;
                          p2 = null;
                        }
                      }
                      find_const = true;
                      continue;
                    }
                  }
                  if ( (A1 <= A2) && (B1 >= B2))
                  // levoe plecho j-oi vklucheno v levoe plecho t-oi
                  {
                    // sluchai 3C
                    if (A2 - A1 > D1 - D2) // naklon vpravo
                    {
                      folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                      int n_down = A2 - A1 + 1;
                      int n_up = C2 - C1 + 1;
                      //chislo par, kotorye nado razrushit
                      if ( (n_up > 0 && n_up >= n_down) || (n_up <= 0))
                      { // razrushaem snizu
                        int n_break = n_down; // chislo par, k-rye nado razrushit'
                        left = helix1.Left.substring(0, n_break);
                        right = helix1.Right.substring(0, n_break);
                        // energiya razrushennogo kusochka
                        energy_break += helix1.energy_calculation(left, right);

                        // energiya ostavshegosya kusochka
                        left = helix1.Left.substring(n_break);
                        right = helix1.Right.substring(n_break);
                        energy_of_part += helix1.energy_calculation(left, right);
                        t_exists[t] = true;
                        pair p = new pair(A2, D2);
                        pairsAD[t] = p;
                        pairsBC[t] = p;
                        // ot j-oi spirali est' kusochek
                        t_exists[j] = true;
                        pair p1 = new pair(A1 + n_break, D1 - n_break);
                        pair p2 = new pair(B1, C1);
                        pairsAD[j] = p1;
                        pairsBC[j] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }
                      if (n_up > 0 && n_up < n_down)
                      {
                        // razrushaem sverhu
                        int n_break = n_up;
                        left = helix1.Left.substring(helix1.Left.length() -
                            n_break);
                        right = helix1.Right.substring(helix1.Left.length() -
                            n_break);
                        energy_break += helix1.energy_calculation(left, right);
                        t_exists[t] = true;
                        pair p = new pair(B2, C2);
                        pairsAD[t] = p;
                        pairsBC[t] = p;
                        // ot j-oi spirali est' kusochek
                        t_exists[j] = true;
                        pair p1 = new pair(A1, D1);
                        pair p2 = new pair(B1 - n_break, C1 + n_break);
                        pairsAD[j] = p1;
                        pairsBC[j] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }

                      // esli nesovmestima tol'ko s odnoi spiral'yu,
                      // to nam nuzhna budet konstanta obratnogo perehoda
                      if (uncompatible == 1)
                      {
                        if (B1 == B2 && D1 <= D2)
                        // znachit svobodnaya iniciaciya nevozmozhna
                        // nuzhno razrushit' odnu paru sverhu
                        {
                          n_back = 1;
                          // energiya razrushennogo kusochka
                          energy_break_back = 0;
                          pair p = new pair(B1, C1);
                          // ot j-oi odna para
                          t_exists_back[j] = true;
                          pairsAD_back[j] = p;
                          pairsBC_back[j] = p;

                          t_exists_back[t] = true;
                          pair p1 = new pair(A2, D2);
                          pair p2 = new pair(B2 - n_back, C2 + n_back);
                          pairsAD_back[t] = p1;
                          pairsBC_back[t] = p2;
                          p = null;
                          p1 = null;
                          p2 = null;
                        }
                        else // kogda vozmozhna svobodnaya iniciaciya
                        {
                          n_down = 0;
                          n_up = 0;
                          pair p = null;
                          if (D2 < D1)
                          {
                            n_down = Math.min( (A2 - A1), (D1 - D2));
                            p = new pair(A1, D1);
                          }
                          if (C1 < C2)
                          {
                            n_up = Math.min( (B1 - B2), (C2 - C1));
                            p = new pair(B1, C1);
                          }
                          t_exists_back[t] = true;
                          pair p1 = new pair(A2, D2);
                          pair p2 = new pair(B2, C2);
                          pairsAD_back[t] = p1;
                          pairsBC_back[t] = p2;
                          n_back = n_down + n_up;
                          if (n_back == 0) System.err.println("n = 0, obr. 3C");
                          // ot j-oi odna para
                          t_exists_back[j] = true;
                          pairsAD_back[j] = p;
                          pairsBC_back[j] = p;
                          p = null;
                          p1 = null;
                          p2 = null;
                        }
                      }
                      find_const = true;
                      continue;
                    }
                    // sluchai 3D
                    if (A2 - A1 < D1 - D2) // naklon vlevo
                    {
                      folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                      //chislo par, kotorye nado razrushit
                      int n_down = D1 - D2 + 1;
                      int n_up = B1 - B2 + 1;
                      if ( (n_down > 0 && n_down >= n_up) || n_down <= 0)
                      {
                        int n_break = n_up; // chislo par, k-rye nado razrushit'
                        left = helix1.Left.substring(helix1.Left.length() -
                            n_break);
                        right = helix1.Right.substring(helix1.Left.length() -
                            n_break);
                        energy_break += helix1.energy_calculation(left, right);

                        // energiya ostavshegosya kusochka
                        left = helix1.Left.substring(0,
                            helix1.Left.length() - n_break);
                        right = helix1.Right.substring(0,
                            helix1.Left.length() - n_break);
                        energy_of_part += helix1.energy_calculation(left, right);

                        t_exists[t] = true;
                        pair p = new pair(B2, C2);
                        pairsAD[t] = p;
                        pairsBC[t] = p;
                        // ot j-oi spirali est' kusochek
                        t_exists[j] = true;
                        pair p1 = new pair(A1, D1);
                        pair p2 = new pair(B1 - n_break, C1 + n_break);
                        pairsAD[j] = p1;
                        pairsBC[j] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }
                      if (n_down > 0 && n_down < n_up)
                      {
                        int n_break = n_down;
                        left = helix1.Left.substring(0, n_break);
                        right = helix1.Right.substring(0, n_break);
                        // energiya razrushennogo kusochka
                        energy_break += helix1.energy_calculation(left, right);

                        // energiya ostavshegosya kusochka
                        left = helix1.Left.substring(n_break);
                        right = helix1.Right.substring(n_break);
                        energy_of_part += helix1.energy_calculation(left, right);
                        t_exists[t] = true;
                        pair p = new pair(A2, D2);
                        pairsAD[t] = p;
                        pairsBC[t] = p;
                        // ot j-oi spirali est' kusochek
                        t_exists[j] = true;
                        pair p1 = new pair(A1 + n_break, D1 - n_break);
                        pair p2 = new pair(B1, C1);
                        pairsAD[j] = p1;
                        pairsBC[j] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }

                      // esli nesovmestima tol'ko s odnoi spiral'yu,
                      // to nam nuzhna budet konstanta obratnogo perehoda
                      if (uncompatible == 1)
                      {
                        if (A1 == A2 && C1 >= C2)
                        // znachit svobodnaya iniciaciya nevozmozhna
                        // nuzhno razrushit' odnu paru snizu
                        {
                          n_back = 1;
                          // energiya razrushennogo kusochka
                          energy_break_back = 0;
                          // ot j-oi odna para
                          pair p = new pair(A1, D1);
                          t_exists_back[j] = true;
                          pairsAD_back[j] = p;
                          pairsBC_back[j] = p;

                          t_exists_back[t] = true;
                          pair p1 = new pair(A2 + n_back, D2 - n_back);
                          pair p2 = new pair(B2, C2);
                          pairsAD_back[t] = p1;
                          pairsBC_back[t] = p2;
                          p = null;
                          p1 = null;
                          p2 = null;
                        }
                        else // kogda vozmozhna svobodnaya iniciaciya
                        {
                          n_down = 0;
                          n_up = 0;
                          pair p = null;
                          if (D2 < D1)
                          {
                            n_down = Math.min( (A2 - A1), (D1 - D2));
                            p = new pair(A1, D1);
                          }
                          if (C1 < C2)
                          {
                            n_up = Math.min( (B1 - B2), (C2 - C1));
                            p = new pair(B1, C1);
                          }
                          t_exists_back[t] = true;
                          pair p1 = new pair(A2, D2);
                          pair p2 = new pair(B2, C2);
                          pairsAD_back[t] = p1;
                          pairsBC_back[t] = p2;
                          n_back = n_down + n_up;
                          if (n_back == 0)
                          {
                            System.err.println("n = 0, obr. 3D");
                          }
                          // ot j-oi odna para
                          t_exists_back[j] = true;
                          pairsAD_back[j] = p;
                          pairsBC_back[j] = p;
                          p = null;
                          p1 = null;
                          p2 = null;
                        }
                      }
                      find_const = true;
                      continue;
                    }
                  }
                  // sluchai 3E
                  // levoe plecho t-oi vklyucheno v pravoe plecho j-oi
                  if (A2 >= C1 && B2 <= D1 && C2 >= D1)
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    //chislo par, kotorye nado razrushit
                    int n_break = D1 - B2 + 1; // chislo par, k-rye nado razrushit'
                    left = helix1.Left.substring(0, n_break);
                    right = helix1.Right.substring(0, n_break);
                    // energiya razrushennogo kusochka
                    energy_break += helix1.energy_calculation(left, right);

                    // energiya ostavshegosya kusochka
                    left = helix1.Left.substring(n_break);
                    right = helix1.Right.substring(n_break);
                    energy_of_part += helix1.energy_calculation(left, right);
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    // ot j-oi spirali est' kusochek
                    t_exists[j] = true;
                    pair p1 = new pair(A1 + n_break, D1 - n_break);
                    pair p2 = new pair(B1, C1);
                    pairsAD[j] = p1;
                    pairsBC[j] = p2;
                    p = null;
                    p1 = null;
                    p2 = null;

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      if (C1 == A2)
                      // znachit svobodnaya iniciaciya nevozmozhna
                      // nuzhno razrushit' odnu paru snizu
                      {
                        n_back = 1;
                        // energiya razrushennogo kusochka
                        energy_break_back = 0;
                        p = new pair(B1, C1);
                        // ot j-oi odna para
                        t_exists_back[j] = true;
                        pairsAD_back[j] = p;
                        pairsBC_back[j] = p;

                        t_exists_back[t] = true;
                        p1 = new pair(A2 + n_back, D2 - n_back);
                        p2 = new pair(B2, C2);
                        pairsAD_back[t] = p1;
                        pairsBC_back[t] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }
                      else // kogda vozmozhna svobodnaya iniciaciya
                      {
                        n_back = A2 - C1;
                        if (n_back == 0)
                        {
                          System.err.println("n = 0, obr. 3E");
                        }

                        t_exists_back[t] = true;
                        p1 = new pair(A2, D2);
                        p2 = new pair(B2, C2);
                        pairsAD_back[t] = p1;
                        pairsBC_back[t] = p2;

                        p = new pair(B1, C1);
                        // ot j-oi odna para
                        t_exists_back[j] = true;
                        pairsAD_back[j] = p;
                        pairsBC_back[j] = p;
                        p = null;
                        p1 = null;
                        p2 = null;

                      }
                    }
                    find_const = true;
                    continue;
                  }
                  // sluchai 3F
                  // pravoe plecho t-oi vklyucheno v levoe plecho j-oi
                  if (C2 >= A1 && D2 <= B1 && B2 <= A1)
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    //chislo par, kotorye nado razrushit
                    int n_break = C2 - A1 + 1; // chislo par, k-rye nado razrushit'
                    left = helix1.Left.substring(0, n_break);
                    right = helix1.Right.substring(0, n_break);
                    // energiya razrushennogo kusochka
                    energy_break += helix1.energy_calculation(left, right);

                    // energiya ostavshegosya kusochka
                    left = helix1.Left.substring(n_break);
                    right = helix1.Right.substring(n_break);
                    energy_of_part += helix1.energy_calculation(left, right);
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    // ot j-oi spirali est' kusochek
                    t_exists[j] = true;
                    pair p1 = new pair(A1 + n_break, D1 - n_break);
                    pair p2 = new pair(B1, C1);
                    pairsAD[j] = p1;
                    pairsBC[j] = p2;
                    p = null;
                    p1 = null;
                    p2 = null;

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      if (B1 == D2)
                      // znachit svobodnaya iniciaciya nevozmozhna
                      // nuzhno razrushit' odnu paru snizu
                      {
                        n_back = 1;
                        energy_break_back = 0;
                        // ot j-oi odna para
                        p = new pair(B1, C1);
                        t_exists_back[j] = true;
                        pairsAD_back[j] = p;
                        pairsBC_back[j] = p;

                        t_exists_back[t] = true;
                        p1 = new pair(A2 + n_back, D2 - n_back);
                        p2 = new pair(B2, C2);
                        pairsAD_back[t] = p1;
                        pairsBC_back[t] = p2;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }
                      else // nichego razrushat' ne nado, svobodnaya iniciaciya
                      {
                        n_back = B1 - D2;
                        t_exists_back[t] = true;
                        p1 = new pair(A2, D2);
                        p2 = new pair(B2, C2);
                        pairsAD_back[t] = p1;
                        pairsBC_back[t] = p2;

                        p = new pair(B1, C1);
                        // ot j-oi odna para
                        t_exists_back[j] = true;
                        pairsAD_back[j] = p;
                        pairsBC_back[j] = p;
                        p = null;
                        p1 = null;
                        p2 = null;
                      }

                    }
                    find_const = true;
                    continue;
                  }
                  //*/
                  //KONETS OBRATNOGO SLUCHAIA 3

                  //===========================================================================
                  // SLUCHAI 4
                  // dlya iniciacii neobhodimo chastichnoe razrushenie
                  // novaya spiral' vytesnit staruyu
                  // 4A - novaia spiral idiot sleva snizu naparavo naverh
                  // schitaya staruiu horisontalnoi, paravie plechi peresekautsia
                  if ( (A2 - A1 < D1 - D2) && // uslovia naklona vlevo
                      (D2 >= C1) && (D2 <= D1) && // peresechenie pravih plechei
                      (C2 <= C1) && B2 <= A1)
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    n += 1;
                    // int n_break = C1 - C2 + 1; // kolichestvo par, kotorye nado razrushit' v nachale spirali
                    int n_break = D1 - D2 + 1; // skol'ko nado razrushit' v nachale spirali
                    left = helix1.Left.substring(0, n_break);
                    right = helix1.Right.substring(0, n_break);
                    // energiya razrushennogo kusochka
                    energy_break += helix1.energy_calculation(left, right);

                    // energiya ostavshegosya kusochka
                    left = helix1.Left.substring(n_break);
                    right = helix1.Right.substring(n_break);
                    energy_of_part += helix1.energy_calculation(left, right);

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(A2, D2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    // ot j-oi spirali est' kusochek
                    t_exists[j] = true;
                    pair p1 = new pair(A1 + n_break, D1 - n_break);
                    pair p2 = new pair(B1, C1);
                    pairsAD[j] = p1;
                    pairsBC[j] = p2;
                    p = null;
                    p1 = null;
                    p2 = null;

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      pair p1_back = null;
                      pair p2_back = null;
                      pair p3_back = null;

                      n_back = 1;
                      n_break = C1 - C2 + 1; // kolichestvo par, kotorye nado razrushit' v konce spirali
                      left = helix2.Left.substring(helix2.Left.length() -
                          n_break);
                      right = helix2.Right.substring(helix2.Left.length() -
                          n_break);
                      energy_break_back = helix2.energy_calculation(left, right);

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      p1_back = new pair(B1, C1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      p2_back = new pair(A2, D2);
                      p3_back = new pair(B2 - n_break, C2 + n_break);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;

                    }

                    find_const = true;
                    continue;
                  }
                  //KONETS 4A

                  // 4B - novaia spiral idet sleva snizu naparavo naverh
                  // schitaya staruiu horisontalnoi, levye plechi peresekautsia
                  if ( (A2 - A1 < D1 - D2) && // uslovia naklona vlevo
                      (B2 >= A1) && (B2 <= B1) && // peresechenie levyh plechei
                      (A2 <= A1) && (D2 <= C1))
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    n += 1;
                    int n_break = B1 - B2 + 1; // kolichestvo par, kotorye nado razrushit' v konce spirali
                    left = helix1.Left.substring(helix1.Left.length() -
                                                 n_break);
                    right = helix1.Right.substring(helix1.Left.length() -
                        n_break);
                    energy_break += helix1.energy_calculation(left, right);

                    // energiya ostavshegosya kusochka
                    left = helix1.Left.substring(0,
                                                 helix1.Left.length() - n_break);
                    right = helix1.Right.substring(0,
                        helix1.Left.length() - n_break);
                    energy_of_part += helix1.energy_calculation(left, right);

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    // ot j-oi spirali est' kusochek
                    t_exists[j] = true;
                    pair p1 = new pair(A1, D1);
                    pair p2 = new pair(B1 - n_break, C1 + n_break);
                    pairsAD[j] = p1;
                    pairsBC[j] = p2;
                    p = null;
                    p1 = null;
                    p2 = null;

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      n_break = A1 - A2 + 1; // kolichestvo par, kotorye nado razrushit' v nachale spirali
                      left = helix2.Left.substring(0, n_break);
                      right = helix2.Right.substring(0, n_break);
                      energy_break_back = helix2.energy_calculation(left, right);

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      pair p1_back = new pair(A1, D1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      pair p2_back = new pair(A2 + n_break, D2 - n_break);
                      pair p3_back = new pair(B2, C2);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;

                    }

                    find_const = true;
                    continue;
                  }
                  //KONETS 4B

                  // 4C - novaia spiral idet sprava snizu nalevo naverh
                  // schitaya staruiu horisontalnoi, levye plechi peresekautsia
                  if ( (A2 - A1 > D1 - D2) && // uslovia naklona napravo
                      (A2 >= A1) && (A2 <= B1) && // peresechenie levyh plechei
                      (B2 >= B1) && (C2 >= D1) /*(C2 > C1)*/)
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    n += 1;
                    int n_break = A2 - A1 + 1; // kolichestvo par, kotorye nado razrushit' v nachale spirali
                    left = helix1.Left.substring(0, n_break);
                    right = helix1.Right.substring(0, n_break);
                    energy_break += helix1.energy_calculation(left, right);

                    // energiya ostavshegosya kusochka
                    left = helix1.Left.substring(n_break);
                    right = helix1.Right.substring(n_break);
                    energy_of_part += helix1.energy_calculation(left, right);

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(A2, D2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    // ot j-oi spirali est' kusochek
                    t_exists[j] = true;
                    pair p1 = new pair(A1 + n_break, D1 - n_break);
                    pair p2 = new pair(B1, C1);
                    pairsAD[j] = p1;
                    pairsBC[j] = p2;
                    p = null;
                    p1 = null;
                    p2 = null;

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      n_break = B2 - B1 + 1; // kolichestvo par, kotorye nado razrushit' v konce spirali
                      left = helix2.Left.substring(helix2.Left.length() -
                          n_break);
                      right = helix2.Right.substring(helix2.Left.length() -
                          n_break);
                      energy_break_back = helix2.energy_calculation(left, right);

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      pair p1_back = new pair(B1, C1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot j-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      pair p2_back = new pair(A2, D2);
                      pair p3_back = new pair(B2 - n_break, C2 + n_break);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;

                    }

                    find_const = true;
                    continue;
                  }
                  //KONETS 4C

                  // 4D - novaia spiral idet sprava snizu nalevo naverh
                  // schitaya staruiu horisontalnoi, pravye plechi peresekautsia
                  if ( (A2 - A1 > D1 - D2) && // uslovia naklona napravo
                      (C2 >= C1) && (C2 <= D1) && // peresechenie pravyh plechei
                      (D2 >= D1) && (A2 >= B1))
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    n += 1;
                    int n_break = C2 - C1 + 1; // kolichestvo par, kotorye nado razrushit' v konce spirali
                    left = helix1.Left.substring(helix1.Left.length() -
                                                 n_break);
                    right = helix1.Right.substring(helix1.Left.length() -
                        n_break);
                    energy_break += helix1.energy_calculation(left, right);

                    // energiya ostavshegosya kusochka
                    left = helix1.Left.substring(0,
                                                 helix1.Left.length() - n_break);
                    right = helix1.Right.substring(0,
                        helix1.Left.length() - n_break);
                    energy_of_part += helix1.energy_calculation(left, right);

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    // ot j-oi spirali est' kusochek
                    t_exists[j] = true;
                    pair p1 = new pair(A1, D1);
                    pair p2 = new pair(B1 - n_break, C1 + n_break);
                    pairsAD[j] = p1;
                    pairsBC[j] = p2;
                    p = null;
                    p1 = null;
                    p2 = null;

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      n_break = D2 - D1 + 1; // skol'ko nado razrushit' v nachale spirali
                      left = helix2.Left.substring(0, n_break);
                      right = helix2.Right.substring(0, n_break);
                      // energiya razrushennogo kusochka
                      energy_break_back = helix2.energy_calculation(left, right);

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      pair p1_back = new pair(A1, D1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      pair p2_back = new pair(A2 + n_break, D2 - n_break);
                      pair p3_back = new pair(B2, C2);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;

                    }

                    find_const = true;
                    continue;
                  }
                  //KONETS 4D

                  // 4E - peresekayutsya raznoimennye plechi (pravoe j-oi i levoe t-oi)
                  if (A2 < C1 && B2 > C1 && B2 <= D1 && C2 > D1)
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    n += 1;
                    int n_break = D1 - B2 + 1; // kolichestvo par, kotorye nado razrushit' v konce spirali
                    left = helix1.Left.substring(0, n_break);
                    right = helix1.Right.substring(0, n_break);
                    energy_break += helix1.energy_calculation(left, right);

                    // energiya ostavshegosya kusochka
                    left = helix1.Left.substring(n_break);
                    right = helix1.Right.substring(n_break);
                    energy_of_part += helix1.energy_calculation(left, right);

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    // ot j-oi spirali est' kusochek
                    t_exists[j] = true;
                    pair p1 = new pair(A1 + n_break, D1 - n_break);
                    pair p2 = new pair(B1, C1);
                    pairsAD[j] = p1;
                    pairsBC[j] = p2;
                    p = null;
                    p1 = null;
                    p2 = null;

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      n_break = C1 - A2 + 1; // kolichestvo par, kotorye nado razrushit' v konce spirali
                      left = helix2.Left.substring(0, n_break);
                      right = helix2.Right.substring(0, n_break);
                      energy_break_back = helix2.energy_calculation(left, right);

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      pair p1_back = new pair(B1, C1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      pair p2_back = new pair(A2 + n_break, D2 - n_break);
                      pair p3_back = new pair(B2, C2);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;
                    }
                    find_const = true;
                    continue;
                  }

                  // 4F - peresekayutsya raznoimennye plechi (levoe j-oi i pravoe t-oi)
                  if (D2 > B1 && C2 < B1 && C2 >= A1 && B2 < A1)
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    n += 1;
                    int n_break = C2 - A1 + 1; // kolichestvo par, kotorye nado razrushit' v konce spirali
                    left = helix1.Left.substring(0, n_break);
                    right = helix1.Right.substring(0, n_break);
                    energy_break += helix1.energy_calculation(left, right);

                    // energiya ostavshegosya kusochka
                    left = helix1.Left.substring(n_break);
                    right = helix1.Right.substring(n_break);
                    energy_of_part += helix1.energy_calculation(left, right);

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    // ot j-oi spirali est' kusochek
                    t_exists[j] = true;
                    pair p1 = new pair(A1 + n_break, D1 - n_break);
                    pair p2 = new pair(B1, C1);
                    pairsAD[j] = p1;
                    pairsBC[j] = p2;
                    p = null;
                    p1 = null;
                    p2 = null;

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      n_break = D2 - B1 + 1; // kolichestvo par, kotorye nado razrushit' v konce spirali
                      left = helix2.Left.substring(0, n_break);
                      right = helix2.Right.substring(0, n_break);
                      energy_break_back = helix2.energy_calculation(left, right);

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      pair p1_back = new pair(B1, C1);
                      pairsAD_back[j] = p1_back;
                      pairsBC_back[j] = p1_back;
                      // ot t-oi spirali est' kusochek
                      t_exists_back[t] = true;
                      pair p2_back = new pair(A2 + n_break, D2 - n_break);
                      pair p3_back = new pair(B2, C2);
                      pairsAD_back[t] = p2_back;
                      pairsBC_back[t] = p3_back;
                      p1_back = null;
                      p2_back = null;
                      p3_back = null;

                    }

                    find_const = true;
                    continue;
                  }

                  // 4H - peresekayutsya oba plecha
                  if ( (A2 - A1 > D1 - D2) && // uslovia naklona napravo
                      (C2 >= C1) && (C2 <= D1) && // peresechenie pravyh plechei
                      (A2 >= A1) && (A2 <= B1)) // peresechenie levyh plechei
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    n += 1;
                    int n_break = 0;
                    int n_up = C2 - C1 + 1;
                    int n_down = A2 - A1 + 1;
                    if (n_up < n_down)
                    {
                      n_break = n_up;
                      left = helix1.Left.substring(helix1.Left.length() -
                          n_break);
                      right = helix1.Right.substring(helix1.Left.length() -
                          n_break);
                      energy_break += helix1.energy_calculation(left, right); // energiya ostavshegosya kusochka

                      // ot t-oi spirali est' tol'ko 1 para
                      t_exists[t] = true;
                      pair p = new pair(B2, C2);
                      pairsAD[t] = p;
                      pairsBC[t] = p;
                      // ot j-oi spirali est' kusochek
                      t_exists[j] = true;
                      pair p1 = new pair(A1, D1);
                      pair p2 = new pair(B1 - n_break, C1 + n_break);
                      pairsAD[j] = p1;
                      pairsBC[j] = p2;
                      p = null;
                      p1 = null;
                      p2 = null;
                    }
                    else
                    {
                      n_break = n_down;
                      left = helix1.Left.substring(0, n_break);
                      right = helix1.Right.substring(0, n_break);

                      energy_break += helix1.energy_calculation(left, right);
                      // ot t-oi spirali est' tol'ko 1 para
                      t_exists[t] = true;
                      pair p = new pair(A2, D2);
                      pairsAD[t] = p;
                      pairsBC[t] = p;
                      // ot j-oi spirali est' kusochek
                      t_exists[j] = true;
                      pair p1 = new pair(A1 + n_break, D1 - n_break);
                      pair p2 = new pair(B1, C1);
                      pairsAD[j] = p1;
                      pairsBC[j] = p2;
                      p = null;
                      p1 = null;
                      p2 = null;
                    }

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      n_up = B2 - B1 + 1;
                      n_down = D2 - D1 + 1;
                      if (n_up < n_down)
                      {
                        n_break = n_up;
                        left = helix2.Left.substring(helix2.Left.length() -
                            n_break);
                        right = helix2.Right.substring(helix2.Left.length() -
                            n_break);
                        energy_break_back = helix2.energy_calculation(left,
                            right); // energiya ostavshegosya kusochka

                        // ot j-oi spirali est' tol'ko 1 para
                        t_exists_back[j] = true;
                        pair p1_back = new pair(B1, C1);
                        pairsAD_back[j] = p1_back;
                        pairsBC_back[j] = p1_back;
                        // ot t-oi spirali est' kusochek
                        t_exists_back[t] = true;
                        pair p2_back = new pair(A2, D2);
                        pair p3_back = new pair(B2 - n_break, C2 + n_break);
                        pairsAD_back[t] = p2_back;
                        pairsBC_back[t] = p3_back;
                        p1_back = null;
                        p2_back = null;
                        p3_back = null;

                      }
                      else
                      {
                        n_break = n_down;
                        left = helix2.Left.substring(0, n_break);
                        right = helix2.Right.substring(0, n_break);
                        energy_break_back = helix2.energy_calculation(left,
                            right); // energiya ostavshegosya kusochka

                        // ot j-oi spirali est' tol'ko 1 para
                        t_exists_back[j] = true;
                        pair p1_back = new pair(A1, D1);
                        pairsAD_back[j] = p1_back;
                        pairsBC_back[j] = p1_back;
                        // ot t-oi spirali est' kusochek
                        t_exists_back[t] = true;
                        pair p2_back = new pair(A2 + n_break, D2 - n_break);
                        pair p3_back = new pair(B2, C2);
                        pairsAD_back[t] = p2_back;
                        pairsBC_back[t] = p3_back;
                        p1_back = null;
                        p2_back = null;
                        p3_back = null;

                      }
                    }
                    find_const = true;
                    continue;
                  }

                  if ( (A2 - A1 < D1 - D2) && // uslovia naklona vlevo
                      (B2 >= A1) && (B2 <= B1) && // peresechenie levyh plechei
                      (D2 >= C1) && (D2 <= D1)) // peresechenie pravih plechei
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    n += 1;
                    int n_break = 0;
                    int n_up = B1 - B2 + 1;
                    int n_down = D1 - D2 + 1;
                    if (n_up < n_down)
                    {
                      n_break = n_up;
                      left = helix1.Left.substring(helix1.Left.length() -
                          n_break);
                      right = helix1.Right.substring(helix1.Left.length() -
                          n_break);
                      energy_break += helix1.energy_calculation(left, right); // energiya ostavshegosya kusochka

                      // ot t-oi spirali est' tol'ko 1 para
                      t_exists[t] = true;
                      pair p = new pair(B2, C2);
                      pairsAD[t] = p;
                      pairsBC[t] = p;
                      // ot j-oi spirali est' kusochek
                      t_exists[j] = true;
                      pair p1 = new pair(A1, D1);
                      pair p2 = new pair(B1 - n_break, C1 + n_break);
                      pairsAD[j] = p1;
                      pairsBC[j] = p2;
                      p = null;
                      p1 = null;
                      p2 = null;
                    }
                    else
                    {
                      n_break = n_down;
                      left = helix1.Left.substring(0, n_break);
                      right = helix1.Right.substring(0, n_break);
                      energy_break += helix1.energy_calculation(left, right);
                      // ot t-oi spirali est' tol'ko 1 para
                      t_exists[t] = true;
                      pair p = new pair(A2, D2);
                      pairsAD[t] = p;
                      pairsBC[t] = p;
                      // ot j-oi spirali est' kusochek
                      t_exists[j] = true;
                      pair p1 = new pair(A1 + n_break, D1 - n_break);
                      pair p2 = new pair(B1, C1);
                      pairsAD[j] = p1;
                      pairsBC[j] = p2;
                      p = null;
                      p1 = null;
                      p2 = null;
                    }

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = 1;
                      n_down = A1 - A2 + 1;
                      n_up = C1 - C2 + 1;
                      if (n_up < n_down)
                      {
                        n_break = n_up;
                        left = helix2.Left.substring(helix2.Left.length() -
                            n_break);
                        right = helix2.Right.substring(helix2.Left.length() -
                            n_break);
                        energy_break_back = helix2.energy_calculation(left,
                            right); // energiya ostavshegosya kusochka

                        // ot j-oi spirali est' tol'ko 1 para
                        t_exists_back[j] = true;
                        pair p1_back = new pair(B1, C1);
                        pairsAD_back[j] = p1_back;
                        pairsBC_back[j] = p1_back;
                        // ot t-oi spirali est' kusochek
                        t_exists_back[t] = true;
                        pair p2_back = new pair(A2, D2);
                        pair p3_back = new pair(B2 - n_break, C2 + n_break);
                        pairsAD_back[t] = p2_back;
                        pairsBC_back[t] = p3_back;
                        p1_back = null;
                        p2_back = null;
                        p3_back = null;
                      }
                      else
                      {
                        n_break = n_down;
                        left = helix2.Left.substring(0, n_break);
                        right = helix2.Right.substring(0, n_break);
                        energy_break_back = helix2.energy_calculation(left,
                            right);
                        // ot j-oi spirali est' tol'ko 1 para
                        t_exists_back[j] = true;
                        pair p1_back = new pair(A1, D1);
                        pairsAD_back[j] = p1_back;
                        pairsBC_back[j] = p1_back;
                        // ot t-oi spirali est' kusochek
                        t_exists_back[t] = true;
                        pair p2_back = new pair(A2 + n_break, D2 - n_break);
                        pair p3_back = new pair(B2, C2);
                        pairsAD_back[t] = p2_back;
                        pairsBC_back[t] = p3_back;
                        p1_back = null;
                        p2_back = null;
                        p3_back = null;
                      }
                    }
                    find_const = true;
                    continue;
                  }

                  //KONETS 4 //*/

                  //===============================================
                  // vozmozhna svobodnaya iniciaciya, no novay spiral' vytesnit staruyu
                  // skol'zyaschuyu petlyu nel'zya obrazovat', potomu chto peresechenie slishkom bol'shoe
                  if ( (A2 - A1 < D1 - D2) && // uslovia naklona vlevo
                      (D1 < D2) && (C2 <= D1) && // pravye plechi peresekayutsya
                      ( (C2 - C1 < rna1.minpartialhelix) || // bez peresecheniya ostaetsya sparennymi men'she (minhelix) par
                       (D2 - D1 < rna1.minpartialhelix)))
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    //chislo par, godnih k initsiatsii spirali snizu ot sush. 1-i (j-oi) spirali
                    n += D2 - D1;

                    //                energy_break += helix1.energy;
                    if (n < 0) System.err.println("n = " + n + "\t3");

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(A2, D2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      if (C2 > C1)
                        n_back = C2 - C1;

                      t_exists_back[t] = true;
                      pair p1 = new pair(A2, D2);
                      pair p2 = new pair(B2, C2);
                      pairsAD_back[t] = p1;
                      pairsBC_back[t] = p2;

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      p = new pair(B1, C1);
                      pairsAD_back[j] = p;
                      pairsBC_back[j] = p;
                      p = null;
                      p1 = null;
                      p2 = null;
                    }

                    find_const = true;
                    continue;
                  }

                  if ( (A2 - A1 < D1 - D2) && // uslovia naklona vlevo
                      (B2 > B1) && (B1 >= A2) && // levye plechi peresekayutsya
                      ( (B1 - B2 < rna1.minpartialhelix) ||
                       (A1 - A2 < rna1.minpartialhelix)))
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    n += B2 - B1;

                    if (n < 0) System.err.println("n = " + n + "\t4");

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      if (A1 < A2)
                      {
                        n_back = A2 - A1;
                      }
                      t_exists_back[t] = true;
                      pair p1 = new pair(A2, D2);
                      pair p2 = new pair(B2, C2);
                      pairsAD_back[t] = p1;
                      pairsBC_back[t] = p2;

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      p = new pair(A1, D1);
                      pairsAD_back[j] = p;
                      pairsBC_back[j] = p;
                      p = null;
                      p1 = null;
                      p2 = null;
                    }

                    find_const = true;
                    continue;
                  }

                  if ( (A2 - A1 > D1 - D2) && // naklon vpravo
                      (A1 > A2) && (B2 >= A1) && // levye plechi peresekayutsya
                      ( (B1 - B2 < rna1.minpartialhelix) ||
                       (A1 - A2 < rna1.minpartialhelix)))
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    //chislo par, godnih k initsiatsii spirali snizu ot sush. 1-i (j-oi) spirali
                    n += A1 - A2;

                    //                 energy_break += helix1.energy;
                    if (n < 0) System.err.println("n = " + n + "\t5");

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(A2, D2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      if (B1 > B2)
                      {
                        n_back = B1 - B2;
                      }
                      t_exists_back[t] = true;
                      pair p1 = new pair(A2, D2);
                      pair p2 = new pair(B2, C2);
                      pairsAD_back[t] = p1;
                      pairsBC_back[t] = p2;

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      p = new pair(B1, C1);
                      pairsAD_back[j] = p;
                      pairsBC_back[j] = p;
                      p = null;
                      p1 = null;
                      p2 = null;
                    }

                    find_const = true;
                    continue;
                  }

                  if ( (A2 - A1 > D1 - D2) && // naklon vpravo
                      (C1 <= D2) && (C2 < C1) && // pravye plechi peresekayutsya
                      ( (C2 - C1 < rna1.minpartialhelix) || // bez peresecheniya ostaetsya sparennymi men'she (minhelix) par
                       (D2 - D1 < rna1.minpartialhelix)))
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    //chislo par, godnih k initsiatsii spirali sverhu ot sush. 1-i (j-oi) spirali
                    n += C1 - C2;

                    //               energy_break += helix1.energy;
                    if (n < 0) System.err.println("n = " + n + "\t6");

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      if (D2 < D1)
                      {
                        n_back = D1 - D2;
                      }
                      t_exists_back[t] = true;
                      pair p1 = new pair(A2, D2);
                      pair p2 = new pair(B2, C2);
                      pairsAD_back[t] = p1;
                      pairsBC_back[t] = p2;

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      p = new pair(A1, D1);
                      pairsAD_back[j] = p;
                      pairsBC_back[j] = p;
                      p = null;
                      p1 = null;
                      p2 = null;
                    }

                    find_const = true;
                    continue;
                  }
                  //==================================================================
                  // spirali idut posledovatel'no (peresekayutsya raznoimennye plechi)
                  if ( (A2 > C1) && (A2 <= D1) && (B2 > D1) && // levoe plecho sp. t peresekaetsya s pravym sp. j
                      ( (A2 - C1 < rna1.minpartialhelix) ||
                       (B2 - D1 < rna1.minpartialhelix)))
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    //chislo par, godnih k initsiatsii
                    n += B2 - D1;

                    //               energy_break += helix1.energy;
                    if (n < 0) System.err.println("n = " + n + "\t6");

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;
                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = A2 - C1;

                      t_exists_back[t] = true;
                      pair p1 = new pair(A2, D2);
                      pair p2 = new pair(B2, C2);
                      pairsAD_back[t] = p1;
                      pairsBC_back[t] = p2;

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      p = new pair(B1, C1);
                      pairsAD_back[j] = p;
                      pairsBC_back[j] = p;

                      p = null;
                      p1 = null;
                      p2 = null;
                    }

                    find_const = true;
                    continue;

                  }
                  if ( (D2 < B1) && (D2 >= A1) && (C2 < A1) && // levoe plecho sp. i peresekaetsya s pravym sp. j
                      ( (B1 - D2 < rna1.minpartialhelix) ||
                       (A1 - C2 < rna1.minpartialhelix)))
                  {
                    folding_next.delete_helix(j); // razrushaem j-uyu spiral'
                    //chislo par, godnih k initsiatsii
                    n += A1 - C2;

                    if (n < 0) System.err.println("n = " + n + "\t6");

                    // ot t-oi spirali est' tol'ko 1 para
                    t_exists[t] = true;
                    pair p = new pair(B2, C2);
                    pairsAD[t] = p;
                    pairsBC[t] = p;
                    p = null;

                    // esli nesovmestima tol'ko s odnoi spiral'yu,
                    // to nam nuzhna budet konstanta obratnogo perehoda
                    if (uncompatible == 1)
                    {
                      n_back = B1 - D2;
                      t_exists_back[t] = true;
                      pair p1 = new pair(A2, D2);
                      pair p2 = new pair(B2, C2);
                      pairsAD_back[t] = p1;
                      pairsBC_back[t] = p2;

                      // ot j-oi spirali est' tol'ko 1 para
                      t_exists_back[j] = true;
                      p = new pair(B1, C1);
                      pairsAD_back[j] = p;
                      pairsBC_back[j] = p;

                      p = null;
                      p1 = null;
                      p2 = null;
                    }

                    find_const = true;
                    continue;

                  }
                  j1 = j;
                  // esli ne podobrali tip nesovmestimosti, to prosto razrushaem
                  folding_next.delete_helix(j);
                }
              if (!find_const || j1 != -1) System.err.println(
                  "Ќе подобрали тип несовместимости " + "t = " + t + " j = " +
                  j1);
//
              // sozdaem perehodnoe sostoyanie
              transit_foldings t_folding = new transit_foldings(pairsAD,
                  pairsBC, t_exists, heliset1);
              // schitatem energiyu perehodnogo sostoyaniya
              transit_energy = t_folding.loops_energy +
                  folding_next.folding_energy - folding_next.loops_energy + // energiya spiralei v strukture s razrushennymi spiralyami
                  energy_of_part; // energiya ne razrushennogo kusochka spirali
              /*              if (sliding_ends > 0) n = heliset1.get_length_without_slid(t,
                                folding_next.exists);
//*/
              folding_next.add_helix(t); // dobavlyaem spiral t
//             energy_break = helix2.energy;

              //             if (n <= 0) System.err.println("n = " + n + "\t7");
              if (n == 0)
              {
                n = 1;
              }
              konst_form = k0 * n * Math.exp( (current_folding.loops_energy -
                                               t_folding.loops_energy +
                                               energy_break) / RT);

              // eto perehod vne gruppy,
              // esli spiral', kotoraya obrazuetsya nesovmestima bolish chem s odnoi spiral'yu
              // ili konstanta men'she poroga
              // i esli ego net vnutri tekuschei gruppy
              if ( (uncompatible > 1) || (konst_form < threshold_form) &&
                  !processing_elements.containsKey(folding_next.exists))
              {
                // zapominaem konstanru perehoda
                current_folding.outside_const[t] = konst_form;
                // i dobavlyaem folding_next v spisok grupp, kak nezapolnennyi element
                // esli ego esche tam net
                if (!rna1.groups_list1.hash.containsKey(folding_next.
                    exists))
                { // my zapominaem etot perehod i ego konstantu
                  current_folding.perehody[t] = folding_next;
                  group gr1 = new group(folding_next, false);
                  rna1.groups_list1.add_new_group(gr1);

                  /*                  rna1.writestruct.Writeln("not completed group");
                                    //rna1.writestruct.Writeln(gr1.toString());
                   rna1.writestruct.Writeln(folding_next.toString(rna1.sequence.
                                        length()));//*/
                }
                else
                {
                  current_folding.perehody[t] = this.find_perehod(folding_next);
                }

                continue;
                //
              }
              // esli nesovmestima tol'ko s odnoi i konstanta bol'she poroga
              // to eto pretendent na chlenstvo v gruppe
              // i nado proveryat' konstantu razrusheniya
              if ( (uncompatible == 1) && (konst_form > threshold_form))
              {
                // ischem konstantu razrusheniyu
                konst_dis = 0;
                // sozdaem obratnoe perehodnoe sostoyanie
                transit_foldings t_folding_back = new transit_foldings(
                    pairsAD_back,
                    pairsBC_back, t_exists_back, heliset1);

                if (t_folding_back == null)
                {
                  System.err.println("t_folding_back == null");
                  return;
                }
                // schitaem obratnuyu konstantu
                konst_dis = k0 * n_back * Math.exp( (folding_next.loops_energy -
                    t_folding_back.loops_energy + energy_break_back) / RT);

                // esli i konstanta razrusheniya bol'she poroga,
                // to eto polnocennyi chlen gruppy i my ego tuda dobavlyaem, esli ego esche ne bylo
                if ( (konst_dis > threshold_dis) &&
                    (!processing_elements.containsKey(folding_next.exists)))
                {

                  // esli takaya structura byla v spiske, to my ee udalyaem
                  if (rna1.groups_list1.hash.containsKey(folding_next.
                      exists))
                  {
                    // esli est', vytaskivaem nomer gruppy
                    Integer group_index = (Integer) rna1.groups_list1.hash.
                        get(folding_next.exists);
                    group gr = (group) rna1.groups_list1.get(group_index.
                        intValue());
                    if (gr != null && gr.completed)
                    {
                      System.err.println("Try to remove completed group (3)" +
                                         group_index.intValue());
                      System.err.println(gr.toString());
                    }
                    // to kladem v elements tot folding, kotoryi byl v toi nezapolnennoi gruppe
                    // t.e. prosto perepisyvaem owner_group elementu gr
                    foldings gr_el = (foldings) gr.processing_elements.get(
                        folding_next.exists);
                    if (gr_el == null) System.err.println("null element");
                    else
                      add_element_to_group(gr_el);

                      // udalyaem etu gruppu iz spiska
                    rna1.groups_list1.remove_group(group_index.intValue());
                    folding_next = null;
                  }
                  else add_element_to_group(folding_next);
                  continue;
                }
                // esli men'she i ego net vnutri gruppy, to tozhe naruzhu
                if (konst_dis < threshold_dis &&
                    !processing_elements.containsKey(folding_next.exists))
                {
                  // zapominaem konstanru perehoda
                  current_folding.outside_const[t] = konst_form;
                  // i dobavlyaem folding_next v spisok grupp, kak nezapolnennyi element
                  // esli ego esche tam net
                  if (!rna1.groups_list1.hash.containsKey(folding_next.
                      exists))
                  { // my zapominaem etot perehod i ego konstantu
                    current_folding.perehody[t] = folding_next;
                    group gr1 = new group(folding_next, false);
                    rna1.groups_list1.add_new_group(gr1);

                    /*                  rna1.writestruct.Writeln("not completed group");
                                      //rna1.writestruct.Writeln(gr1.toString());
                     rna1.writestruct.Writeln(folding_next.toString(rna1.sequence.
                                          length()));//*/
                  }
                  else
                  {
                    current_folding.perehody[t] = this.find_perehod(
                        folding_next);
                  }
                  continue; //
                }
                t_folding_back = null;
              }
              t_folding = null;
              folding_next = null;
              pairsAD = null;
              pairsAD_back = null;
              pairsBC = null;
              pairsBC_back = null;
              t_exists_back = null;
              t_exists = null;
//              folding_next.up_loop = null;
//              folding_next.down_loop = null;
            }
          }
        }
      }
      else // esli spiral' mozhet razrushit'sya
      { // t.e. folding.exists[t]=true i eta spiral mozhet razrushit'sya

        double[] get_energy_length = heliset1.
            get_length_energi_without_slid(t, current_folding.exists);
        double split_energy = get_energy_length[1];
        int n = (int) get_energy_length[0];

        folding_next.delete_helix(t); // razrushaem spiral i
        // konstanta razrusheniya
        konst_dis = (n - 1) * k0 * Math.exp(split_energy / (RT));
        // schitaem konstantu skorosti obrazovaniya (zavisit tol'ko ot energii petel' i dliny spirali!!!!)
        konst_form = k0 * (n - 1) *
            Math.exp( (folding_next.loops_energy -
                       current_folding.loops_energy) / (RT));
        // sravnivaem s porogom obe konstanty, esli oni obe bol'she poroga
        // to eto perehod vnutri gruppy
        if ( (konst_dis > threshold_dis) && (konst_form > threshold_form))
        {
          // esli ego tam esche net
          if (!processing_elements.containsKey(folding_next.exists))
          {
            // esli takaya structura byla v spiske, to my ee udalyaem
            if (rna1.groups_list1.hash.containsKey(folding_next.exists))
            {
              // esli est', vytaskivaem nomer gruppy
              Integer group_index = (Integer) rna1.groups_list1.hash.get(
                  folding_next.exists);
              group gr = (group) rna1.groups_list1.get(group_index.
                  intValue());
              if (gr != null && gr.completed)
              {
                System.err.println("Try to remove completed group (4)" +
                                   group_index.intValue());
                System.err.println(gr.toString());
                System.err.println("Current group");
                System.err.println(this.toString());
              }
              // to kladem v elements tot folding, kotoryi byl v toi nezapolnennoi gruppe
              // t.e. prosto perepisyvaem owner_group elementu gr
              foldings gr_el = (foldings) gr.processing_elements.get(
                  folding_next.exists);
              if (gr_el == null) System.err.println("null element");
              else
                add_element_to_group(gr_el);
                // udalyaem etu gruppu iz spiska
              rna1.groups_list1.remove_group(group_index.intValue());
              folding_next = null;
            }
            else add_element_to_group(folding_next);
            continue;
          }
        }
        /*       else // esli hotya by odna men'she, to eto perehod vne gruppy
               {
                 // i my zapominaem etot perehod i ego konstantu
                 gr_el.outside_const[t] = konst_dis;
                 gr_el.perehody[t] = folding_next;
                 // i dobavlyaem v spisok grupp, kak ne zapolnennyi element
                 // esli ego esche tam net i esli ego net vnutri tekuschei gruppy
         if (!structure.groups_list1.hash.containsKey(folding_next.exists) &&
                     (!processing_elements.containsKey(folding_next.exists)))
                 {
                   group_element gr_el2 = new group_element(folding_next);
                   group gr1 = new group(gr_el2);
                   structure.groups_list1.add_new_group(gr1);

         /*              rna1.writestruct.Writeln("not completed group");
                                   //rna1.writestruct.Writeln(gr1.toString());
           rna1.writestruct.Writeln(folding_next.toString(rna1.sequence.
                                       length()));
                   }
                   continue;
                 }//*/

      }
//      folding_next.up_loop = null;
//      folding_next.down_loop = null;
    }
    //  System.err.println("GC in collect_group");
//    System.gc();

    // govorim chto etot folding prosmotren
    current_folding.checked = true;
  }

  public String toString()
  {
    StringBuffer sb = new StringBuffer();
    if (elements != null)
      for (int i = 0; i < elements.length; i++)
      {
        sb.append("element #" + i + "\n" +
                  elements[i].toString(rna1.sequence.length()));
      }
    else
    {
      if (processing_elements != null)
      {
        Set key_set = processing_elements.keySet(); // vytaskivaem set klyuchei
        Iterator iter;
        for (iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
        {
          boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
          // i prisvaevaem ocherednome elementu massiva
          foldings gr_el = (foldings) processing_elements.get(
              key_exists);
          sb.append(gr_el.toString(rna1.sequence.length()));
        }
      }
    }
    sb.append("completed = " + this.completed + "\n");
    return sb.toString();
  }

  public String XML(int id)
 {
   StringBuffer sb = new StringBuffer();
   sb.append("<Group ID = \"" + id + "\" Completed = \"" + completed + "\" BestElementID = \"" + best + "\">" + "\n");
   if (elements != null)
     for (int i = 0; i < elements.length; i++)
     {
  //     sb.append( elements[i].XML(rna1.sequence.length(), i) + "\n");
     }
   sb.append("</Group>" + "\n");
   return sb.toString();
 }


  //=====================================
  public String toString(int length)
  {
    StringBuffer sb = new StringBuffer();
    sb.append("time = " + (live_time/rna1.M) + "\n" + elements[best].toString(length));
    if (rna1.print_graph == true)
    {
      sb.append("graph 1: ");
      // delim vse elementy massiva graph na M, chtoby poluchit' veroyatnost'
      for (int i = 0; i < rna1.graph_step; i++)
      {
        // delim na kolichestvo prohodov
        graph_1[i] = (graph_1[i] / (rna1.M));
        // delim na shag, chtoby poluchilas' veroyatnost'
        graph_1[i] = (graph_1[i] / (rna1.Teta/rna1.graph_step));
        sb.append(graph_1[i] + ", ");
      }
    }
    sb.append("\n" + print_helises(elements[best]));
    return sb.toString();
  }

  public String print_helises(foldings f)
  {
    StringBuffer sb = new StringBuffer();
    Vector true_elem = f.index_exists_true();
    // pechataem vybrannye elementy, t.e. suschestvuyuschie spirali
    for (int j = 0; j < true_elem.size(); j++)
    {
      Integer J = (Integer) true_elem.get(j);
      sb.append(J + "\n");
    }
   // if (true_elem.isEmpty()) sb.append("\n");
    return sb.toString();
  }



  //=======================================
  // nahodit element gruppy po exists i vydaet ego nomer v elements
  int find_element(boolean[] exists)
  {
    if (elements != null)
      for (int i = 0; i < elements.length; i++)
      {
        Comparator_exists com = new Comparator_exists();
        if (com.compare(elements[i].exists, exists) == 0)
          return i;
      }
    return -1;
  }

//===================================================
  // zapolnyaem grafik dlya structury
  void add_to_graph(double CurrTime, double LiveTime)
  {
    // kladem vremya v grafiki dly structury (graph 1) i spiralei (graph 2)
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
          // grafik struktury
        graph_1[rna1.graph_step - 1] += (delta);
        break;
      }
      if (next_time > time_end) // esli konec zhizni lezhit v predydushem elemente
      {
        // sootvetstvuyuschii otrezok vremeni dobavlyaem k elementu i
        // k grafiku struktury
        graph_1[i] += (time_end - time_start);
      }
      else // esli konec zhizni lezhit dal'she
      {
        // sootvetstvuyuschii otrezok vremeni dobavlyaem k elementu i
        // v grafike struktury
        graph_1[i] += (next_time - time_start);
        // i sdvigaem vremya nachala
        time_start = next_time;
      }
    }
  }
}


  /*//=================================================
  // vydaet massiv konstant perehoda do momenta, naskol'ko vyrosla cep' (po indeksu massiva exists)
  double[] get_consts(int exist_index)
  {
    double[] curr_consts = null;
    boolean curr_key[] = new boolean[heliset_size]; // maksimal'nyi klyuch
    // formiruem maksimal'nyi klyuch
    // vse, chto men'she tekuschei pozicii = true,
    // vse, chto bol'she = false
    for (int i = 0; i<= exist_index; i++) curr_key[i] = true;
    for (int i = exist_index+1; i < heliset_size; i++) curr_key[i] = false;
    // elementy s kluchom men'shim curr_key
    TreeMap curr_elements = (TreeMap) elements.headMap(curr_key);
    // kolichestvo elementov, s kluchom men'shim curr_key
     int curr_size = curr_elements.size();
    // zadaem razmer massiva konstant, kak (exist_index*curr_size)
    curr_consts = new double[((exist_index+1)*curr_size)];

    Set key_set = elements.keySet(); // vytaskivaem set klyuchei
    int k = 0; // nomer
    for (Iterator iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
    {
      boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
      group_element gr_el = (group_element) curr_elements.get(key_exists); // izvlekaem sootvetstvuyuschii element
      for (int i = 0; i <= exist_index; i++)
        curr_consts[(exist_index+1)*k + i] = gr_el.outside_const[i];
      k++;
    }
    return curr_consts;
  }//*/





//==============================================================================================================
// hranit gruppy structur, ob'edinennye vmeste, esli konstanty perehodov men'she group_threshold
// kazhdyi element predstavlyaet soboi ob'ekt classa group
public class groups_list extends Vector
{
//  Vector groups; // vektor iz ob'ektov classa group
  TreeMap hash; // po exists govorit k kakoi gruppe otnositsya (nomer v vektore gpoups)

  public groups_list()
  {
    Comparator_exists C_exists; // klyuch, po kotoromu budut uporyadocheny elementy hash
    C_exists = new Comparator_exists();
    hash = new TreeMap(C_exists);
  }

  void add_new_group(group gr)
  {
    // dobavili ee
    this.add(gr);
    // nashli nomer etoi gruppy v spiske
    Integer current_index = new Integer(this.size() - 1);
    if (gr.elements != null)
      // dlya vseh elementov gruppy zapisyvaem nomer gruppy v hash
      // i samim foldingam, nomer gruppy
      for (int i = 0; i < gr.elements.length; i++)
      {
        hash.put(gr.elements[i].exists, current_index);
      }
    else // esli elements pustoe, znachit gruppa ne zapolnena

    // i v hash kladem processing_elements
    if (gr.processing_elements != null)
    {
      Set key_set = gr.processing_elements.keySet(); // vytaskivaem set klyuchei
      Iterator iter;
      for (iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
      {
        boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
        // i kladem v hash
        foldings gr_el = (foldings) gr.processing_elements.get(key_exists);
        hash.put(gr_el.exists, current_index);
      }
    }
    gr.group_index = size()-1;
  }


  // dobavit' v konkretnot mesto
  void add_new_group(int index, group gr)
  {
    // udalyaem s etogo mesta element
    this.remove(index);
    // dobavili ee na mesto index
    this.add(index, gr);
    if (gr != null)
    if (gr.elements != null)

      // dlya vseh elementov gruppy zapisyvaem nomer gruppy v hash
      // i samim foldingam, nomer gruppy
      for (int i = 0; i < gr.elements.length; i++)
      {
        hash.put(gr.elements[i].exists, new Integer(index));
      }
    else // esli elements pustoe, znachit gruppa ne zapolnena

    // i v hash kladem processing_elements
    if (gr.processing_elements != null)
    {
      Set key_set = gr.processing_elements.keySet(); // vytaskivaem set klyuchei
      Iterator iter;
      for (iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
      {
        boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
        // i kladem v hash
        foldings gr_el = (foldings) gr.processing_elements.get(key_exists);
        hash.put(gr_el.exists, new Integer(index));
      }
    }
  if (gr == null)
    System.err.println("Try to add empty group " + index);
    gr.group_index = index;
  }


  void remove_group(int index)
  {
    // snachala udalim iz hash
    group gr = (group) get(index);
    if (gr != null)
    // esli gruppa zapolnena
    // to iz hash udalyaem elements
      if (gr.completed)
        for (int i = 0; i < gr.elements.length; i++)
        {
          hash.remove(gr.elements[i].exists);
        }
      else // esli net, to udalyaem processing_elements
      if (gr.processing_elements != null)
      {
        Set key_set = gr.processing_elements.keySet(); // vytaskivaem set klyuchei
        Iterator iter;
        for (iter = key_set.iterator(); iter.hasNext(); ) // i perebiraem ih iteraterom
        {
          boolean[] key_exists = (boolean[]) iter.next(); // vytaskivaem ocherednoi klyuch'
          // i kladem v hash
          foldings gr_el = (foldings) gr.processing_elements.get(
              key_exists);
          hash.remove(gr_el.exists);
        }
      }
    if (gr == null)
    System.err.println("Try to remove empty group " + index);
    // potom i samu gruppu
    this.remove(index);
    // a na ee mesto kladem pustoi ob'ekt
    this.add(index, null);
  }

//
  public String toString(boolean complete)
  {
    StringBuffer sb = new StringBuffer();
    if (this.size() != 0)
    {
      sb.append("SIZE = " + this.size()+ "\n");
      for (int i = 0; i < size(); i++)
      {
        group gr = (group) get(i);
        if (gr != null)
          if (complete && gr.completed)// pechatat' tol'ko zapolnennye gruppy
            sb.append("GROUP # = " + i + "\n" + gr.toString());
        if (!complete)// pechatat' vse gruppy
          sb.append("GROUP # = " + i + "\n" + gr.toString());
      }
    }
    return sb.toString();
  }


  public String XML(boolean complete)
 {
   StringBuffer sb = new StringBuffer();
   sb.append("<GroupList>");
   sb.append("<SIZE>" + this.size()+ "</SIZE>" + "\n");
   if (this.size() != 0)
   {
     for (int i = 0; i < size(); i++)
     {
       group gr = (group) get(i);
       if (gr != null)
         if (complete && gr.completed)// pechatat' tol'ko zapolnennye gruppy
           sb.append(gr.XML(i) + "\n");
       if (!complete)// pechatat' vse gruppy
         sb.append(gr.XML(i) + "\n");
     }
   }
   sb.append("</GroupList>");
   return sb.toString();
 }


}
