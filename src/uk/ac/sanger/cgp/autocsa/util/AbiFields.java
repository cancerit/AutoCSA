package uk.ac.sanger.cgp.autocsa.util;

public class AbiFields {
  public static String [] labels={"DATA", "DATA", "DATA", "DATA","DATA",
                                  "DATA", "DATA", "DATA", "DATA","DATA",
                                  "DyeN", "DyeN", "DyeN", "DyeN","DyeN",
                                  "MODL"};
  public static int [] indices={1,  2,  3,  4, 105,
                                9, 10, 11, 12, 205,
                                1,  2,  3,  4,   5,
                                1};
  public AbiFields() {
  }
}
/*
//--- define data to read from abi sample file via Tags
lable = ['DATA', 'DATA', 'DATA', 'DATA','DATA',$ ; rawn data
        'DATA', 'DATA', 'DATA', 'DATA','DATA',$ ; analysed data
        'SPAC', 'DYEP',            $ ; spacing & mob file name
        'DyeN', 'DyeN', 'DyeN', 'DyeN',$ ; four dye names
        'FWO_', 'S/N%',            $ ; base order & sig strength
        'CMNT', 'SMPL',            $ ; sample name &
        'MCHN', 'LANE',            $ ; mach name & capillary number
        'PBAS', 'PLOC'             $ ; called bases & their locations
       ]

//--- define corresponding abi sample file Tag indices
list = [1, 2, 3, 4,105,     $      ; list for raw data
       9, 10, 11, 12,205,   $      ; list for analysed data
       1, 1,         $      ;
       1, 2, 3, 4,   $      ; list for four different dyes
       1, 1,         $
       1, 1,         $
       1, 1,         $
       2, 1          $
       ]
//--- End data definition
*/
