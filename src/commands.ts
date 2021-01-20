import { command } from "commander";

interface ribxz_command {
  flag         : string;
  longform_flag: string;
  description  : string;
  arguments?   : Array<string>;
}

type Commands = 
'inspect' |
'overview'

export const availableCommands: Record<Commands, ribxz_command> = {
  inspect: {
    description:
      "Parse ribosomal file and inspect properties (according to additional flags)",
    flag: "-insp",
    longform_flag: "--inspect",
    arguments: ['pdbid', 'protnum']
  },
    //    commander apparently implements a `help` method on the main module. Silly. Calling this overview.
   overview: {
    description: "Display help for the application",
    flag: "-overview",
    longform_flag: "--overview",
  },
};


// want my own cl-parser kinda badly

// export const get_activations = (
//   availcommands: Record<string, ribxz_command>
// ) => {
//   var activations: Record<string, boolean> = {};
//   Object.entries(AvailableCommands)
//     .map(command => {
//       return { [command[0]]: false };
//     })
//     .map(activation => {
//       activations = { ...activations, ...activation };
//     });

//   return activations;
// };
