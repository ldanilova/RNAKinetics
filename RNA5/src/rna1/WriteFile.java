package rna1;

/**
 * <p>Title: ReadFile</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Luda
 * @version 1.0
 */

//$Id: WriteFile.java,v 1.10 2003/10/29 01:35:55 favorov Exp $

public class WriteFile
{
  private String Name;
  public WriteFile(String FileName, boolean if_append)
  {
    Name=FileName;
    try
    {
      java.io.FileWriter wdata = new java.io.FileWriter(Name,if_append);
      // f=true means append;
      // f=false - write anew
    }
    catch (Exception e)
    {
      System.err.println("Cannot open the output file.");
    }
  }
  // write one string

  public void Writeln(String s)
  {
    try
    {
      java.io.FileWriter wdata = new java.io.FileWriter(Name,true);
      // f=true means append;
      // f=false - write anew
      wdata.write(s);
      wdata.write("\n");
      wdata.close();
    }
    catch (Exception e)
    {
      System.err.println("File is not create");
    }
  }

  public void Write(String s)
 {
   try
   {
     java.io.FileWriter wdata = new java.io.FileWriter(Name,true);
     // f=true means append;
     // f=false - write anew
     wdata.write(s);
     wdata.close();
   }
   catch (Exception e)
   {
     System.err.println("File is not create");
   }
 }

  // write array of strings
  public void Write(Object o)
  {
    Writeln(o.toString());
  }

  public void Write(Object[] vo)
  {
    for (int i=0;i<vo.length;i++)
      Writeln(vo.toString());
  }

  public void Write(java.util.Vector v)
  {
    for (int i=0;i<v.size();i++)
      Writeln(i + "\n"+ v.get(i).toString());
  }

  public void Write(double d)
  {
    Double D = new Double(d);
    Writeln(D.toString());
  }

  public void Close(java.io.FileWriter output)
  {
    try
    {
      output.close();
    }
    catch (Exception e)
    {
      System.err.println("Cannot close the output file.");
    }
  }
}
