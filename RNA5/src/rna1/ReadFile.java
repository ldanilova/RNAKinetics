package rna1;

import java.util.*;
import java.io.*;

/**
 * <p>Title: ReadFile</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Luda
 * @version 1.0
 */

//$Id: ReadFile.java,v 1.6 2003/10/27 19:16:26 favorov Exp $

public class ReadFile
{
  Vector sequences = new Vector();
  Vector names = new Vector();
  Vector structures = new Vector();


  public ReadFile(String Filename)
  {
   // System.out.println("Hi!");
    String curr_line = "";

    File indata = new File(Filename);
    try
    {
      BufferedReader bdata = new BufferedReader(new FileReader(indata));
      try
      {
        while (bdata.ready())
        {
          curr_line = bdata.readLine(); // prochitali stroku iz faila
          if (curr_line.equals("")) continue; // esli pustaya
          int tab_counter = 0;
          String name = "", seq = "", struct = "";
          for (int i = 0; i<curr_line.length(); i++)
          {
            String curr_char = curr_line.substring(i, i+1);
            if (tab_counter==0)
            { // t.e. chitaem name
              if (!curr_char.equals("\t")) name = name.concat(curr_char);
              else
                {
                  tab_counter++;
                  continue;
                }
            }
            if ( tab_counter==1 )
            {// chitaem posledovatel'nost'
              if (!curr_char.equals("\t")) seq = seq.concat(curr_char);
              else
                {
                  tab_counter++;
                  continue;
                }
            }
            if ( tab_counter==2 )
           {// chitaem strukturu
             if (!curr_char.equals("\t")) struct = struct.concat(curr_char);
             else
               {
                 tab_counter++;
                 continue;
               }
           }
          }
          sequences.add(seq);
          names.add(name);
          structures.add(struct);
        }
        bdata.close();
        if (sequences.size()!= names.size() || names.size()!= structures.size())
          System.err.println("Something wrong with data!");
      }
      catch (Exception e)
      {
       System.err.println("No data!");
       WriteFile writelog = new WriteFile(rna1.ID + ".log", false);
       writelog.Writeln("No data!");
      }
    }
    catch (FileNotFoundException e)
    {
      System.err.println("File is not exist");
      WriteFile writelog = new WriteFile(rna1.ID + ".log", false);
      writelog.Writeln("File is not exist");
      System.exit(1);
    }
  }

  public String[] get_sequences()
  {
    String[] seq = new String[sequences.size()];
    sequences.toArray(seq);
    return seq;
  }
  public String[] get_names()
  {
    String[] nam = new String[names.size()];
    names.toArray(nam);
    return nam;
  }
  public String[] get_structures()
  {
    String[] str = new String[structures.size()];
    structures.toArray(str);
    return str;
  }

}
