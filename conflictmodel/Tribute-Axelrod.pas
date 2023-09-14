program Tribute40;
{Ver0.1 begun 10/21/92, Ver 2.0 begun 12/21/92, Ver 3.0 begun 1/13/93}
{Documemtion in Act 66}
{v3.0 Loyalty}
{v3.1 Contiguous alliances, eliminate range of tribute}
{           Event_Output for report_person }
{v3.2 bombs}
{v3.3 Loyalty periodic output, built on v3.1old}
{       2/11 after run 23: move calc_alliance from conduct_fight to make_response to correct w_def}
{v3.4 bombs}
{v3.5 add B5, toughness, }
{        2/12 fix div 0 bug in A2 and A3}
{v3.6 add contiguity=false option to avoid requiring contiguity of alliances}
{v3.7 add loyalty_option for who may increase loyalty: 0=all, 1=either rare, 2=neither rare }
{v3.8 add A4 rule: demand of random other (for checking loyalty only variant)}
{       add zero_init_loyalty: T=Loy(ij)=0 initally as before, F=random initial loyalties}
{      add B6 pay or fight at random}
{v4.0 add tribute counts in final_output and trib_output}
{     add more initial loyalty options }


    const
        Version = 4.0;                                      {Program Version Number}
        test = false;                   {if test=true, no run_number  is read or written.}
                                        {Run # is set to 0 for test=true runs}
        Common_rule_A = 2;                      {rare rule occurs once, common rule elsewhere: A is demander}
                    {   A1=att weakest alliance,  A2=max value*vulnerabilty, A3 is mod A2}
                    {   A4 demand of random (reachable) other}
        Common_rule_B = 4;                      {B is defender's response}
                    {   B1=never pay, B2=always pay, B3=pay if cheaper for self, B4, fight if cheaper for self}
                    {   B5= pay if  fight_cost > toughness*payment}
        Rare_rule_A = 2;                            {Rare rules used only for non-adaptive pops}
        Rare_rule_B = 4;
        Year_max = 100;                             {number of years in a run, often 100 or 1000}
        pop_max = 10;                               {number of populations in a run, usually 50 when adapt=false}
        Report_Person = 5;                  {Usually 5. Report events when person=i,j,aider or helper; no report if  = 0}
        old_random_seed = 0;                {0 means new seed, else enter an old seed to reuse}
        periodic_report_freq = 25;                  {0 means never, 1 means each year, 10 means every 10th yr, etc.}
        contiguity = true;                          {T=contiguity required of alliances, F= contiguity not required}
        loyalty_option = 0;                         {who may increase loyalty: 0=all, 1=either rare, 2=neither rare}
        init_loyalty = 3;                               {0:loyalty_pc(i,j)=0 initially as before, 1:random; for i<>j}
{                                                               2: 0 but L (4,5) = 10%, 3: 0 but L(3,5)=10%}
{Rarely changed constants:     }
        triangle_initial_wealth = false;                {T=old method, F= W_base + or - 100}
        Common_rule_C = 1;                      {C not used: retained for later expansion}
        Common_rule_D = 1;                      {D not used: retained for later expansion}
        Rare_rule_C = 1;                            {C not used: retained for later expansion}
        Rare_rule_D = 1;                            {D not used: retained for later expansion}
        rare_indiv = 5;                             {placement of rare indiv in the pop, eg at i=5. }
        adapt = false;                                  {T=each pop adapts from prev pop, F= pops stay same}
        Run_Number_File_Name = 'Run Number File';
        Imax = 10;                                      {number of actors/pop; 10 for experiment}
        Demand_phase_max = 3;                   {number of possible demands / year; usually 3 when imax=10}
        W_base = 400.0;                             {basic initial wealth, often 400.0}
        Standard_Demand = 250.0;                {demand, often 250.0}
        Destructiveness = 0.25;                 {% other's wealth lost in a fight, often 0.25}
        Allele_Max = 6;                             {no. of possible rules(alleles) of each type(gene)= max of max_allele}
        Gene_Max = 2;                               {number of genes, ie A,B,gives 2}
        mutatable_gene_max = 2;                 {0 means all, n means only first n genes can mutate: in Set_Rare_Rule}
        toughness = 1.5;                                {used in B5. If 1.0 then B5=B3.}
        loyalty_pc_increment = 10;              {change in % loyalty with tribute or fight, should be divisible by 100}
        report_rebellion = false;                   {whether to report potential rebellions in Conduct_Fight}
    type
        actor_type = 1..imax;
        gene_type = 1..4;                               {4 genes: A, B, C , D - but C and D are inactive, need for chrom}
        actor_array = array[actor_type] of real;
        bool_actor_array = array[actor_type] of boolean;
        distance_array = array[actor_type] of 0..imax;
        chrom_type = array[actor_type, gene_type] of 1..Allele_Max;
        offset_type = -imax..imax;                      {offset of target from i}
        max_allele_type = array[1..9] of integer;   {number of possible alleles of a given gene}
        loyalty_pc_type = array[actor_type, actor_type] of 0..100; {Loyalty of i to j in percentag}
        integer_maxtrix_type = array[actor_type, actor_type] of integer;  {for tribute counts}

    var
        initial_datetime, end_datetime: datetimerec;    {date, etc.}
        run_number: integer;
        datafile: text;                                     {input data}
        i, j: actor_type;                                   {index of actors}
        Productivity: actor_array;                  {actors' productiivty}
        W: actor_array;                                 {actor's Wealth, current}
        W_initial: actor_array;                     {   initial wealth}
        W_if_isolated: actor_array;                 {    wealth if isolated, current, ie only productivity changes}
        Demand: real;                                       {demand being made by i on j}
        demand_phase: integer;                          {current demand phase}
        year: integer;                                      {current year}
        income: real;                                       {income from production}
        decide_to_demand: boolean;                  {False=no demand, True=i Demands from j}
        Agree_To_Pay: boolean;                          {False=no pay, True= j will pay i}
        tribute: real;                                      {Tribute paid by j to i}
        final_wealth: text;                             {Output of final wealth file, for Excel}
        periodic_wealth: text;                          {Output of periodic wealth file, for Excel}
        periodic_loyalty: text;                         {Output of periodic loyalty file, for Excel}
        adapt_wealth: text;                             {Output of adaptation wealth file, for Excel, if adapt=true}
        invaders_wealth: text;                          {Output of invaders' wealth file, for Excel, if adapt=true}
        event_wealth: text;                             {Output of events involving the report_person}
        periodic_trib: text;                            {Output of count of per. tribute count matrix}
        random_seed: integer;
        Chrom: chrom_type;                              {e.g. chrom[i,1]=3 means i uses A3}
        rare: bool_actor_array;                     {T=actor is rare type, F=actor is common type}
        offset_target: offset_type;                     {location of target relative to i}
        pop: integer;                                       {current pop number}
        w_att: real;                                        {wealth available for an  attack}
        w_def: real;                                        {Wealth available for a defense}
        pop_wealth: real;                                   {total current wealth of the population}
        fights_since_last_report: integer;          {num fights since the last periodic report}
        years_since_last_report: integer;           {num years since last periodic report}
        max_allele: max_allele_type;                                {number of alleles of a given gene}
        common_rule_label, rare_rule_label: integer;    {for output}
        common_indiv: actor_type;                   {loc of a common indiv, set to be rare_indiv+1}
        start_time, end_time, duration: longint;    {for calc of run's duration}
        initial_hour, end_hour: longint;
        invaders_so_far: integer;                       {numer of invaders so far,when adapt=true}
        duration_since_invasion: integer;           {number of pops since last invasion}
        loyalty_pc: loyalty_pc_type;                    {loyalty of i to j, percentage between 0 and 100}
        aid_att: bool_actor_array;                      {T if k aids attacker}
        help_def: bool_actor_array;                 {T if k helps defender}
        j_temp_reachable: boolean;                  {T if this j_temp would be reachable from i}
        trib_count: integer_maxtrix_type;           {# payments from row to col }
        trib_periodic_count: integer_maxtrix_type; {# payments from row to col in current period}
  { ---------------------------------------------------------------  }
    function min (x, y: real): real;
    begin
        if x < y then
            min := x
        else
            min := y;
    end;
  { ---------------------------------------------------------------  }
    function min_integer (x, y: integer): integer;
    begin
        if x < y then
            min_integer := x
        else
            min_integer := y;
    end;
  { ---------------------------------------------------------------  }
    function random_one_to_n (n: longint): longint;
                            {proc returns a random number between 1 and n; modified from Bennett's random_range}
        var
            ub, lb: integer;
            r: integer;
    begin                                                       {random gives # betw -32768 and 32767}
        ub := 32767 - (32767 mod n);
        lb := -32768 - (-32768 mod n);                  {truncate distrib on 2 ends so that later mod is OK}
        repeat
            r := random;
        until (r <= ub) and (r >= lb);                      {make sure random genrated is in truncated (even) distrib}
        random_one_to_n := abs(r mod n) + 1;
    end;            {random function}

 { ---------------------------------------------------------------  }
    procedure Report_Periodic_Output;
        var
            i, j: actor_type;
    begin
        pop_wealth := 0;
        for i := 1 to imax do
            pop_wealth := pop_wealth + W[i];
        write(periodic_wealth, pop : 3, '   ', year : 3, '  ', fights_since_last_report : 5, '  ', pop_wealth : 6 : 1, '    ');
        write(periodic_wealth, w[1] : 6 : 1, '  ', w[2] : 6 : 1, '  ', w[3] : 6 : 1, '  ', w[4] : 6 : 1, '  ', w[5] : 6 : 1, '  ', w[6] : 6 : 1, '  ');
        writeln(periodic_wealth, w[7] : 6 : 1, '    ', w[8] : 6 : 1, '  ', w[9] : 6 : 1, '  ', w[10] : 6 : 1);
        for i := 1 to imax do                               {Loyalty matrix first}
            begin{do row as outer loop}
                write(periodic_loyalty, pop : 4, '  ', year : 4, '  ', w[i] : 6 : 1, '  ', i : 4, ' ');
                for j := 1 to imax do
                    begin
                        write(periodic_loyalty, loyalty_pc[i, j] : 4, ' ');
                        if j = imax then
                            writeln(periodic_loyalty);
                    end;{j}
                if i = imax then
                    writeln(periodic_loyalty);
            end;{i}
        fights_since_last_report := 0;              {reset for  interval to next periodic report}
        years_since_last_report := 0;               {ditto}
        if year > 0 then                                            {Triubte maxtix next}
            begin
                for i := 1 to imax do                   {do row as outer loop}
                    begin
                        write(periodic_trib, pop : 4, ' ', year : 4, '  ', w[i] : 6 : 1, '  ', i : 4, ' ');
                        for j := 1 to imax do
                            begin
                                write(periodic_trib, trib_periodic_count[i, j] : 4, '   '); {i,j is ok here}
                                trib_periodic_count[i, j] := 0;   {immeidate re-initialize}
                                if j = imax then
                                    writeln(periodic_trib);
                            end;{j}
                        if i = imax then
                            writeln(periodic_trib);
                    end;{i}
            end;{year>0}
    end;
 { ---------------------------------------------------------------  }
    procedure make_person_report;
{Report whenever Report_person indiv is role A, B, aider of attacker, helper of defender}
{From Make_Payment and Conduct_Fight}
        type
            report_string_type = packed array[1..imax] of char;
        var
            report_string: report_string_type;
            k: actor_type;
    begin
        if Agree_to_Pay then
            begin                                                       {if tribute was paid}
                for k := 1 to imax do
                    begin
                        report_string[k] := '-';                                {set to . if nothing applies}
                        if k = i then
                            report_string[k] := 'R';                            {R for receiver of tribute}
                        if k = j then
                            report_string[k] := 'P';                            {P for payor of tribute}
                    end;
            end{if tribute was paid}
        else
            begin
                for k := 1 to imax do                               {if fight}
                    begin
                        report_string[k] := '-';                                {set to . if nothing applies}
                        if aid_att[k] then
                            report_string[k] := 'a';                                {a for aid attacker}
                        if help_def[k] then
                            report_string[k] := 'd';                                {d for help defender}
                        if k = i then
                            report_string[k] := 'A';                            {A for attacker}
                        if k = j then
                            report_string[k] := 'D';                            {D for defender}
                    end;{fight}
            end;{if Agree_to_Pay}
        write(event_wealth, pop : 3, '  ', year : 3, '  ', i : 3, ' ', j : 3, ' ', report_string : 12, '    ');
        write(event_wealth, w[1] : 6 : 1, ' ', w[2] : 6 : 1, '  ', w[3] : 6 : 1, '  ', w[4] : 6 : 1, '  ', w[5] : 6 : 1, '  ', w[6] : 6 : 1, '  ');
        writeln(event_wealth, w[7] : 6 : 1, '   ', w[8] : 6 : 1, '  ', w[9] : 6 : 1, '  ', w[10] : 6 : 1);
    end;
 { ---------------------------------------------------------------  }
    procedure Write_Adapt_and_Invaders_Headers;
    begin
        Write(adapt_wealth, ' Run ', Run_number : 4, '.  Ver ', Version : 4 : 2, ' of ', initial_datetime.month : 2, '/', initial_datetime.day : 2);
        Writeln(adapt_wealth, ' at ', initial_datetime.hour : 2, ':', initial_datetime.minute : 2, '  ');
        Writeln(adapt_wealth, common_rule_label : 5, '  Initial common rule: A, B, C, D');
        Writeln(adapt_wealth, toughness : 7 : 2, '  toughness');
        Writeln(adapt_wealth, pop_max : 3, '    pops');
        Writeln(adapt_wealth, year_max : 3, '   years / pop ');
        Writeln(adapt_wealth, imax : 3, '   indivs / pop ');
        Writeln(adapt_wealth, demand_phase_max : 3, '   demands/year');
        Writeln(adapt_wealth, W_base : 7 : 2, ' initial wealth');
        Writeln(adapt_wealth, standard_demand : 7 : 2, '    standard demand');
        Writeln(adapt_wealth, destructiveness : 7 : 2, '    destructiveness');
        Writeln(adapt_wealth, loyalty_pc_increment : 7, '   loyalty % increment');
        if mutatable_gene_max > 0 then
            writeln(adapt_wealth, mutatable_gene_max : 3, ' initial  genes eligible for mutation')
        else
            writeln(adapt_wealth, ' All genes eligible for mutation');
        writeln(adapt_wealth);
        Writeln(adapt_wealth, 'Col A:     Population number.');
        Writeln(adapt_wealth, 'Col B-E , Common Rule ');
        Writeln(adapt_wealth, 'Col F-I: Rare Rule.');
        Writeln(adapt_wealth, 'Col J, K: Ave Common Wealth, Rare'' s Wealth ');
        Writeln(adapt_wealth, 'Col L:  1 if rare invades, 0 if not');
        Writeln(adapt_wealth, 'Col M, N : Ave Common Range , Rare'' s Range');
        Writeln(adapt_wealth, 'Col O, P: % Common that did better than if isolated, same for rare.');
        writeln(adapt_wealth);
        Writeln(adapt_wealth, 'pop  cA  cB  cC  cD  rA  rB  rC  rD  AvComW  RareW   Inv?    C(L+R)  R(L+R)  C%>isoW R>isoW');

        Write(invaders_wealth, ' Run ', Run_number : 4, '.  Ver ', Version : 4 : 2, ' of ', initial_datetime.month : 2, '/', initial_datetime.day : 2);
        Writeln(invaders_wealth, ' at ', initial_datetime.hour : 2, ':', initial_datetime.minute : 2, '  ');
        Writeln(invaders_wealth, common_rule_label : 5, '   Initial common rule: A, B, C, D');
        Writeln(invaders_wealth, toughness : 7 : 2, '   toughness');
        Writeln(invaders_wealth, pop_max : 3, ' pops');
        Writeln(invaders_wealth, year_max : 3, '    years / pop ');
        Writeln(invaders_wealth, imax : 3, '    indivs / pop ');
        Writeln(invaders_wealth, demand_phase_max : 3, '    demands/year');
        Writeln(invaders_wealth, W_base : 7 : 2, '  initial wealth');
        Writeln(invaders_wealth, standard_demand : 7 : 2, ' standard demand');
        Writeln(invaders_wealth, destructiveness : 7 : 2, ' destructiveness');
        Writeln(invaders_wealth, loyalty_pc_increment : 7, '    loyalty  % increment');
        if mutatable_gene_max > 0 then
            writeln(invaders_wealth, mutatable_gene_max : 3, '  initial  genes eligible for mutation')
        else
            writeln(invaders_wealth, '  All genes eligible for mutation');
        writeln(invaders_wealth);
        writeln(invaders_wealth, 'Col A:      Invaders so far.');
        writeln(invaders_wealth, 'Col B:      Duration, i.e. pops since last invasion.');
        Writeln(invaders_wealth, 'Col C:     Population number.');
        Writeln(invaders_wealth, 'Col D-G , Common Rule ');
        Writeln(invaders_wealth, 'Col H-K: Rare Rule.');
        Writeln(invaders_wealth, 'Col L, M: Ave Common Wealth, Rare'' s Wealth ');
        Writeln(invaders_wealth, 'Col N, O: Ave Common Range , Rare'' s Range');
        Writeln(invaders_wealth, 'Col P, Q: % Common that did better than if isolated, same for rare.');
        writeln(invaders_wealth);
        Writeln(invaders_wealth, 'inv   dur pop cA  cB  cC  cD  rA  rB  rC  rD  AvComW  RareW   C(L+R)  R(L+R)  C%>isoW R>isoW');
    end;{Write_Adapt_and_Invaders_Headers}
  { ---------------------------------------------------------------  }
    procedure Initialize_run;
        var
            i: actor_type;
    begin
        max_allele[1] := 3;                         {number of possible alleles of gene 1, ie A}
        max_allele[2] := 5;                         {number of possible alleles of gene 2, ie B}
        max_allele[3] := 1;                         {number of possible alleles of gene 3, ie C}
        max_allele[4] := 1;                         {number of possible alleles of gene 4, ie D}
        if rare_indiv >= imax then
            begin
                writeln('Warning. Fatal error: rare_indiv may not be last indiv, ie imax');
                halt;
            end;
        common_indiv := rare_indiv + 1;                 {the location of a common indiv, for copying its chrom}
        rewrite(final_wealth, 'Final_Output');          {open output of final wealth file for writing to Excel}
        rewrite(periodic_wealth, 'Periodic_Output');    {open output of periodic wealth file for writing to Excel}
        rewrite(periodic_loyalty, 'Loyalty_Output');    {open output of periodic loyalty file for writing to Excel}
        rewrite(periodic_trib, 'Tribute_Output');       {open output of periodic tribute file for writing to Excel}
        if adapt = true then
            begin
                invaders_so_far := 0;
                duration_since_invasion := 0;
                rewrite(adapt_wealth, 'Adapt_Output'); {open output of adapt wealth file for writing to Excel}
                rewrite(invaders_wealth, 'Invader_Output'); {open output of invaders wealth file for writing to Excel}
            end;
        if report_person <> 0 then
            begin
                rewrite(event_wealth, 'Event_Output');  {open output of event wealth for file for writing to Excel}
{XXXX write header here?}
            end;{report_person<>0}
        gettime(initial_datetime);
        initial_hour := initial_datetime.hour;              {to force long int}
        start_time := 60 * 60 * initial_hour + 60 * initial_datetime.minute + initial_datetime.second;
        Writeln('Output of Axelrod''s Tribute Program, Version ', Version : 5 : 2);
        Write('   This run begun on ', initial_datetime.month : 2, '/', initial_datetime.day : 2);
        Write('/', initial_datetime.year : 4, ' at ', initial_datetime.hour : 2, ':', initial_datetime.minute : 2, '.');
        if old_random_seed = 0 then
            begin                                                   {generate new seed}
                random_seed := initial_datetime.hour + initial_datetime.minute + initial_datetime.second + (initial_datetime.second * 300) + (initial_datetime.minute * initial_datetime.hour) + (initial_datetime.minute * initial_datetime.second);
                randseed := random_seed;
                Writeln(' New random seed ', randseed : 6, '.');
                Writeln(final_wealth, 'Tribute Program. Final Output.', ' New random seed', randseed : 8, '.');
                if periodic_report_freq <> 0 then
                    begin
                        Writeln(periodic_wealth, 'Tribute Program, Periodic Output. ', ' New random seed', randseed : 8, '.');
                        Writeln(periodic_loyalty, 'Tribute Program, Loyalty Output. ', ' New random seed', randseed : 8, '.');
                        Writeln(periodic_trib, 'Tribute Program, Tribute Count Output. ', ' New random seed', randseed : 8, '.');
                    end;
                if report_person <> 0 then
                    Writeln(event_wealth, 'Tribute Program, Event Output. ', ' New random seed', randseed : 8, '.');
                if adapt then
                    begin
                        Writeln(adapt_wealth, 'Tribute Program, Adapt Output. ', ' New random seed', randseed : 8, '.');
                        Writeln(invaders_wealth, 'Tribute Program, Invaders Output. ', ' New random seed', randseed : 8, '.');
                    end;
            end
        else
            begin                                                   {use old seed, which was inputed as constant}
                randseed := old_random_seed;
                Writeln(' Old random seed ', randseed : 6, '.');
                Writeln(final_wealth, 'Tribute Program. Final Output.', ' Old random seed', randseed : 8, '.');
                if periodic_report_freq <> 0 then
                    begin
                        Writeln(periodic_wealth, 'Tribute Program, Periodic Output. ', ' Old random seed', randseed : 8, '.');
                        Writeln(periodic_loyalty, 'Tribute Program, Loyalty Output. ', ' Old random seed', randseed : 8, '.');
                        Writeln(periodic_trib, 'Tribute Program, Tribute Count Output. ', ' Old random seed', randseed : 8, '.');
                    end;
                if report_person <> 0 then
                    Writeln(event_wealth, 'Tribute Program, Event Output. ', ' Old random seed', randseed : 8, '.');
                if adapt then
                    begin
                        Writeln(adapt_wealth, 'Tribute Program, Adapt Output. ', ' Old random seed', randseed : 8, '.');
                        Writeln(invaders_wealth, 'Tribute Program, Invaders Output. ', ' Old random seed', randseed : 8, '.');
                    end;
            end;
        if test then                                             {run number}
            run_number := 0
        else
            begin                   {not a test;  read run #}
                reset(datafile, Run_Number_File_Name);
                case (ioresult) of
                    -43, 17, 19, 21, 24: 
                        begin
                            writeln('File  error opening the run_number file.  ');
                            writeln('This is a fatal error.  Check file name and path in input file and try again ');
                            halt;
                        end;
                    otherwise
                        begin
                            writeln('Run number file opened OK.');
                        end;
                end;                {case}
                readln(datafile, run_number);
                run_number := run_number + 1;
                close(datafile);
                rewrite(datafile, Run_Number_File_Name);
                writeln(datafile, run_number);
                writeln('Run number', run_number : 4, '.');
                close(datafile);
            end;    {if for run number }
        Common_rule_label := 1000 * common_rule_A + 100 * common_rule_B + 10 * common_rule_C + Common_rule_D;
        Rare_rule_label := 1000 * rare_rule_A + 100 * rare_rule_B + 10 * rare_rule_C + rare_rule_D;
        Write(final_wealth, ' Run ', Run_number : 4, ':  Ver ', Version : 4 : 2, ' of ', initial_datetime.month : 2, '/', initial_datetime.day : 2);
        Writeln(final_wealth, ' at ', initial_datetime.hour : 2, ':', initial_datetime.minute : 2, '  ');
        Write(event_wealth, ' Run ', Run_number : 4, ':  Ver ', Version : 4 : 2, ' of ', initial_datetime.month : 2, '/', initial_datetime.day : 2);
        Writeln(event_wealth, ' at ', initial_datetime.hour : 2, ':', initial_datetime.minute : 2, '  ');
        Writeln(event_wealth, '   Key to roles: R=receiver, P=payer of tribute.  A=attacker, a=aider of attacker. D=defender, d=helper of def');
        Writeln(final_wealth, common_rule_label : 5, '  Initial common rule: A, B, C,  D');
        if adapt then
            Writeln(final_wealth, '  Adaptation run, so no initial rare rule.')
        else
            Writeln(final_wealth, rare_rule_label : 5, '    rare rule: A, B, C, D');
        Writeln(final_wealth, toughness : 7 : 2, '  toughness');
        Writeln(final_wealth, pop_max : 3, '    pops');
        Writeln(final_wealth, year_max : 3, '   years / pop ');
        Writeln(final_wealth, imax : 3, '   indivs / pop ');
        Writeln(final_wealth, demand_phase_max : 3, '   demands/year');
        Writeln(final_wealth, W_base : 7 : 2, ' initial wealth');
        Writeln(final_wealth, standard_demand : 7 : 2, '    standard demand');
        Writeln(final_wealth, destructiveness : 7 : 2, '    destructiveness');
        Writeln(final_wealth, loyalty_pc_increment : 7, '   loyalty % increment');
        writeln(final_wealth, contiguity : 7, ' Contiguity of alliances required: T/F');
        writeln(final_wealth, loyalty_option : 3, ' Loyalty dynamics: 0=all, 1=either rare, 2=neither rare');
        writeln(final_wealth, init_loyalty : 3, '   init_loyalty. 0 is init L(i,j)=0. 1 is random; for i<>j. 2 is 0 but L(4,5)=10%. 3 is 0 but L(3,5)=10%.');
            { No blank lines left for future growth}
        writeln(final_wealth, ' pop rare    i   Final Wealth    W if Isolated');
        writeln(event_wealth, 'pop  year    A   B   roles   W1      W2      W3      W4      W5      W6      W7      W8      W9      W10');
        if report_rebellion then
            writeln('Cases with loyalty(i,j)= 100 BEFORE fight, and Wi 0 then
            begin
                Write(periodic_wealth, ' Run ', Run_number : 4, ': Ver ', Version : 4 : 2, ' of ', initial_datetime.month : 2, '/', initial_datetime.day : 2);
                Writeln(periodic_wealth, ' at ', initial_datetime.hour : 2, ':', initial_datetime.minute : 2, '  ');
                writeln(periodic_wealth, 'pop   year    fights  Wpop    W1      W2      W3      W4      W5      W6      W7      W8      W9      W10');
                Write(periodic_loyalty, ' Run ', Run_number : 4, ': Ver ', Version : 4 : 2, ' of ', initial_datetime.month : 2, '/', initial_datetime.day : 2);
                Writeln(periodic_loyalty, ' at ', initial_datetime.hour : 2, ':', initial_datetime.minute : 2, '  ');
                writeln(periodic_loyalty, 'Loyalty, in percent. From row to col.');
                writeln(periodic_loyalty, 'pop  year    W       row 1   2   3   4   5   6   7   8   9   10');
                Write(periodic_trib, ' Run ', Run_number : 4, ': Ver ', Version : 4 : 2, ' of ', initial_datetime.month : 2, '/', initial_datetime.day : 2);
                Writeln(periodic_trib, ' at ', initial_datetime.hour : 2, ':', initial_datetime.minute : 2, '  ');
                writeln(periodic_trib, 'Tribute  Count: No. Times Row Paid Tribute to Col. in Each Period.');
                writeln(periodic_trib, 'pop year    W   row 1   2   3   4   5   6   7   8   9   10');
            end;
        for i := 1 to imax do               {initialize chormosome}
            begin
                if (i = rare_indiv) and (adapt = false) then
                    begin                                           {one rare rule at i=rare_indiv, only if not adapt}
                        chrom[i, 1] := rare_rule_A;
                        chrom[i, 2] := rare_rule_B;
                        chrom[i, 3] := rare_rule_C;
                        chrom[i, 4] := rare_rule_D;
                        rare[i] := True;
                    end
                else                                                    {others are the common rule}
                    begin
                        chrom[i, 1] := common_rule_A;
                        chrom[i, 2] := common_rule_B;
                        chrom[i, 3] := common_rule_C;
                        chrom[i, 4] := common_rule_D;
                        rare[i] := False;
                    end;{if i=rare_indiv, ie rare}
            end;{i loop}
        rare[rare_indiv] := true;                                       {set i=rare_indiv rare, even if adapt, since it will be later}
        if adapt = true then
            Write_Adapt_and_Invaders_Headers;
    end; {initialize run procedure}

  { ---------------------------------------------------------------  }
    procedure Set_Rare_Rule;                            {from initialize_pop for adaptive runs}
        var                                                     {   does mutation.}
            gene_temp: gene_type;                           {  Note: selection is in Conduct_and_Report_Adaptation}
            mutated_gene: gene_type;
            mutated_allele: integer;
            eligible_gene: gene_type;                       {number of genes eligible for mutation}
    begin
        for gene_temp := 1 to gene_max do
            begin
                chrom[rare_indiv, gene_temp] := chrom[common_indiv, gene_temp]; {rare rule starts like a common rule}
            end;{gene_temp}
        if mutatable_gene_max > 0 then              {control variable, if 0 all genes are eligible}
            eligible_gene := mutatable_gene_max
        else
            eligible_gene := gene_max;
        mutated_gene := random_one_to_n(eligible_gene);
        mutated_allele := random_one_to_n(max_allele[mutated_gene] - 1); {-1 to avoid no change}
        if (mutated_allele >= chrom[rare_indiv, mutated_gene]) then {adjust what's being mutated}
            mutated_allele := mutated_allele + 1;                   {   to avoid no change}
        chrom[rare_indiv, mutated_gene] := mutated_allele;
{writeln('Test Set_Rare_Rule, pop, mutated gene, mutated allele', pop : 4, mutated_gene : 4, mutated_allele);}
    end;
  { ---------------------------------------------------------------  }
    procedure Initialize_pop;
        var
            i, j: actor_type;
    begin
        fights_since_last_report := 0;                  {for periodic output under fights column, counted in conduct_fight}
        years_since_last_report := 0;
        pop_wealth := 0;
        for i := 1 to imax do
            begin
                Productivity[i] := 20.0;                        {equal producitivities (gain/artor/year) }
                if triangle_initial_wealth then             {control constant}
                    W[i] := 100 * i                                 {triangular initial wealth distribution}
                else
                    W[i] := W_base + random / 327.67;   {flat initial wealth distribution + or - uniform 100}
                W_initial[i] := W[i];                           {Save initial wealths}
                W_if_isolated[i] := W[i];                       {Start calc of wealth if isolated}
                pop_wealth := pop_wealth + W[i];            {cumulate pop wealth}
                if (i = rare_indiv) and adapt then
                    Set_Rare_Rule;                              {new rare rule set if adaptive run, ie mutation}
                for j := 1 to imax do
                    begin                                           {initialize loyalties}
                        if i = j then
                            loyalty_pc[i, j] := 100                         {loyal to self}
                        else {i<>J}
                            begin
                                if init_loyalty = 0 then                    {if init_loyalty = 0 then loyal to no one else}
                                    loyalty_pc[i, j] := 0;
                                if init_loyalty = 1 then                    {if init_loyalty = 1 then random, eg 0, 10,20..100%  }
                                    begin
                                        loyalty_pc[i, j] := loyalty_pc_increment * (random_one_to_n(1 + 100 div loyalty_pc_increment) - 1);
                                    end; {init_loyalty = 1}
                                if init_loyalty = 2 then                    {if init_loyalty=2 then 0 except L(4,5)=10%) }
                                    begin
                                        if ((i = 4) and (j = 5)) or ((i = 5) and (j = 4)) then
                                            loyalty_pc[i, j] := 10
                                        else
                                            loyalty_pc[i, j] := 0;
                                    end; {initi loyalty = 2}
                                if init_loyalty = 3 then                    {if init_loyalty=3 then 0 except L(3,5)=10%) }
                                    begin
                                        if ((i = 3) and (j = 5)) or ((i = 5) and (j = 3)) then
                                            loyalty_pc[i, j] := 10
                                        else
                                            loyalty_pc[i, j] := 0;
                                    end; {initi loyalty = 3}
                            end;{i<>j}
                    end; {j}
            end; {i loop}
        for i := 1 to imax do
            begin
                for j := 1 to imax do
                    begin                                           {initialize tribute count}
                        trib_count[i, j] := 0;
                        trib_periodic_count[i, j] := 0;
                    end; {j}
            end; {i loop}
        if periodic_report_freq <> 0 then
            begin
                year := 0;
                Report_Periodic_Output;
            end;
    end; {initialize_pop procedure}
  { ---------------------------------------------------------------  }
    procedure Select_Active_Actor;
    begin
        i := random_one_to_n(imax);
    end; {procedure}
  { ---------------------------------------------------------------  }
    procedure Calc_Alliance (var attacker, defender: actor_type; offset: integer);
        label
            2, 3;                                               {offset is from attacker to defender, -is left, + is right}
        var
            other: actor_type;                              {potential alliance member}
            att_al_cand_dis: actor_type;                {distance from i}
            def_cand_dis: actor_type;
            other_offset: offset_type;
            dir_at_all: -1..1;                              {direction for expanding attacking alliance, -1=L, 1=R}
            att_al_cand_offset: integer;                    {from i to attacker alliance candidate}
            dir_to_target: integer;                     {-1 is left, 1 is right}
    begin
        w_att := w[attacker];                       {start with own strength}
        w_def := w[defender];
        for other := 1 to imax do                   {initialize to no help or aid}
            begin
                aid_att[other] := False;
                help_def[other] := False;
            end;
        aid_att[attacker] := True;              {attacker is part of attack}
        help_def[defender] := True;         {defender is part of def}
        if contiguity = true then
            begin                                   {old method: require contiguity of alliances}
                j_temp_reachable := True;           {assume j_temp is reachable until proven otherwise}
                if offset < 0 then
                    dir_to_target := -1                 {def is left of att}
                else
                    dir_to_target := 1;
{Example: ,i = 7, Stage = 1( look left ) , cand_distance = 4 }
{so j_temp = 7 - 4 = 3.}
{if year = 49 then writeln('  Calc_Alliance test: defender=', defender : 3, '. offset=', offset : 3);}
        {XXX test}
                if abs(offset) > 1 then         {there is some space between attacker and j_temp, ie defender}
                    begin
                        for att_al_cand_dis := 1 to abs(offset) - 1 do {check 6,5,4 in att}
                            begin
                                att_al_cand_offset := dir_to_target * att_al_cand_dis;
                                other := (attacker + att_al_cand_offset + imax - 1) mod imax + 1;
                                if loyalty_pc[other, attacker] > loyalty_pc[other, defender] then
                                    begin                                   {contribute to attcker}
                                        w_att := w_att + loyalty_pc[other, attacker] * w[other] / 100;
                                        aid_att[other] := True;                 {other aids attacker in this alliance situtation}
{if year = 49 then writeln ( '    Calc_Alliance test: add ' , other : 3 , ' to attacker alliance on near side.' );}
{XXX}
                                    end
                                else                      {this other won't join attacker, so need go no further in this direction}
                                    begin
                                        j_temp_reachable := false;
                                        goto 3;
{reject this "other" since att all isn 't contig to it}
            {but would still need to check other j_temps further out in this direction}
                                    end;{else}
                            end;{for att_al_cand_dis}
                    end; {abs ({offset)>1, ie some space between attacker and defender}
                for def_cand_dis := 1 to imax - 1 do            {check 2,1...in def}
                    begin
                        other_offset := dir_to_target * def_cand_dis;
                        other := (defender + other_offset + imax - 1) mod imax + 1;
                        if loyalty_pc[other, attacker] < loyalty_pc[other, defender] then
                            begin                                   {contribute to defender}
                                w_def := w_def + loyalty_pc[other, defender] * w[other] / 100;
                                help_def[other] := True;                    {other does  help defender}
{if year = 49 then writeln('    Calc_Alliance test: add ', other : 3, ' to def alliance.');}
{XXX}
                            end
                        else
                            goto 2; { stop search for def }
                    end;{def_cand_dis}
2:
                dir_at_all := -dir_to_target;       {check 8,9... in att}
                for att_al_cand_dis := 1 to imax - 1 do
                    begin
                        att_al_cand_offset := dir_at_all * att_al_cand_dis;
                        other := (attacker + att_al_cand_offset + imax - 1) mod imax + 1;
                        if loyalty_pc[other, attacker] > loyalty_pc[other, defender] then
                            begin
                                w_att := w_att + loyalty_pc[other, attacker] * w[other] / 100;
                                aid_att[other] := True;                 {other aids attacker in this alliance situtation}
{if year = 49 then writeln('    Calc_Alliance test: add ', other : 3, 'to attacker alliance on far side.');}
{XXX}
                            end
                        else
                            goto 3; { stop search in this dir for att }
                    end;{ for att_al_cand_dis}
3:
                ;
            end {if contiguity = true}
        else
            begin                                       { contiguity = false}
                for other := 1 to imax do
                    begin  {already initialized F, except att aids self , def helps self}
                        if (other <> attacker) and (other <> defender) then
                            begin                               {other isn't attacker or defender}
                                if loyalty_pc[other, attacker] > loyalty_pc[other, defender] then
                                    begin
                                        w_att := w_att + loyalty_pc[other, attacker] * w[other] / 100;
                                        aid_att[other] := True;
                                    end; {add attacker}
                                if loyalty_pc[other, defender] > loyalty_pc[other, attacker] then
                                    begin
                                        w_def := w_def + loyalty_pc[other, defender] * w[other] / 100;
                                        help_def[other] := True;
                                    end; {help def}
                            end; {other isn't attakcker or defender}
                    end; {other loop}
            end; {contiguity = false}
    end;{Calc_Alliance}
{-----------------------------------------------------------}
    procedure Calc_Demand_Value (j_temp: actor_type; var value: real);
{called by Make_Demand}
        var
            expect_to_pay: 0..1;                            {0 if not, 1 if so. Used in A3}
            vulnerability: real;
            payment: real;
            targets_upper_cost: real;                   {what target would suffer it rich enough}
            targets_cost_if_fight: real;                    {what  target will actually suffer if there's a fight}
    begin
        case chrom[i, 1] of             {Rules just determine what counts as "value" of a j_temp}
            1:                                  {RULE A1: select weakest alliance in range, ie max difference of W's}
                begin
                    value := W_att - W_def; {value is difference of strength of alliances}
{writeln('    A1: value = ', value : 7 : 2);}
{XXX}
                end;{case A1}
            2:                                  {RULE A2: max payment*vulnerablity}
                begin
                    payment := min(standard_demand, W[j_temp]); {expected payment}
                    if W_att = 0 then                   {prevent div 0}
                        vulnerability := -1.0               {neg vul leads to neg value, so not chosen}
                    else
                        vulnerability := (W_att - W_def) / W_att;           {relative strength}
                    value := payment * vulnerability;
                end;{case A2}
            3:                                  {RULE A3: max value*vulnerability, subject to other's}
                                            {               having it cheaper to pay than fight     }
                begin
                    payment := min(standard_demand, W[j_temp]); {expected payment}
                    if W_att = 0 then                   {prevent div 0}
                        vulnerability := -1.0               {neg vul leads to neg value, so not chosen}
                    else
                        vulnerability := (W_att - W_def) / W_att;           {relative strength}
                    if w_def = 0 then
                        targets_upper_cost := 0         {to avoid div by 0}
                    else
                        targets_upper_cost := destructiveness * W_att * (W[j_temp] / W_def);
                    targets_cost_if_fight := min(targets_upper_cost, W[j_temp]);
                    if payment < targets_cost_if_fight then         {if cheaper to pay}
                        expect_to_pay := 1
                    else
                        expect_to_pay := 0;
                    value := payment * vulnerability * expect_to_pay;
{if (year = 39) and (i = 5) and (j_temp = 4) then}
{begin}
{writeln('    A3: j_temp=', j_temp : 3, 'payment=', payment : 7 : 2, '.  vuln = ', vulnerability : 7 : 2, ' . value = ', value : 7 : 2);}
{writeln('       Wjtemp=', W[j_temp] : 7 : 2, 'W_def = ', W_def : 7 : 2);}
{writeln('       targets_upper_cost', targets_upper_cost : 7 : 2);}
{end;}
{XXX}
                end; {case A3}
            4:                              {RULE A4: Random demand: Give each (reachable) other a positive random value.  }
                begin                       {               The one with highest value will be target.There will always be a target.}
                    value := random + 32768 {gives number from 0 to 64k}
                end;{case of A4}
        end;{case of Make_Demand_Rule}
    end;{Calc_Demand_Value}
  { ---------------------------------------------------------------  }
    procedure Make_Demand;
        label
            1;
        var
            max_value: real;
            value: real;
            j_temp: actor_type;                         {candidate (temp) target}
            j_temp_offset: integer;                                 {- is to left, + is to right of attack}
            targets_upper_cost: real;                   {what target would suffer it rich enough}
            targets_cost_if_fight: real;                    {what  target will actually suffer if there's a fight}
            stage: integer;                                 {1= look left, 2=look right}
            cand_distance: 1..imax;                     {abs distance i to j_temp}
            direction_to_target: integer;               {-1 is left, 1 is right}
    begin
        decide_to_demand := False;                      {no demand if nothing good found}
        max_value := 0;                                     {dont demand if no attractive target}
        if contiguity = true then
            begin
                for stage := 1 to 2 do  {1=left, 2=right}
                    begin
                        j_temp_offset := 0;
                        direction_to_target := 2 * stage - 3;   {1 gives -1, 2 gives +1}
                        for cand_distance := 1 to imax - 1 do   {can go all the way around the circle}
                            begin
                                j_temp_offset := j_temp_offset + direction_to_target;   {go one further away from i}
                                j_temp := (i + j_temp_offset + imax - 1) mod imax + 1;
{if year = 49 then writeln('Make_Demand: j_temp_offset:', j_temp_offset : 3, '. j_temp', j_temp : 3);}
{XXX}
                                calc_alliance(i, j_temp, j_temp_offset); {get reachability, aid,help given this i, j_temp, direction_to_target, cand_offset}
{if year = 49 then writeln('   Make_Demand: reachable:', j_temp_reachable);}
{XXX}
                                if j_temp_reachable then
                                    begin
                                        Calc_Demand_Value(j_temp, value);
                                        if (value > max_value) then         {current j_temp is best so far}
                                            begin
                                                max_value := value;         {reset best seen}
                                                j := j_temp;
                                                offset_target := j_temp_offset;
                                                decide_to_demand := True;                           {decide to demand since value>0}
{if year = 49 then write(' Make_Demand , best value so far : ', value : 6 : 2);}
{XXX}
{if year = 49 then writeln('W_att,W_def', W_att : 7 : 2, ' ', W_def : 7 : 2);}
{XXX}
                                            end; {if value>max value}
                                    end;{if j_temp_reachable}
                                if loyalty_pc[j_temp, i] = 0 then
                                    goto 1 {search no further in this direction}
                            end;{for cand_distance}
1:
                        ;
                    end;{stage}
            end {contiguity = true}
        else
            begin                                               {contiguity = false}
                for j_temp := 1 to imax do
                    begin
                        if j_temp <> i then
                            begin
                                j_temp_offset := 0;     {not used at all here since contiguity = false}
                                calc_alliance(i, j_temp, j_temp_offset);
                                Calc_Demand_Value(j_temp, value);
                                if (value > max_value) then         {current j_temp is best so far}
                                    begin
                                        max_value := value;         {reset best seen}
                                        j := j_temp;
                                        offset_target := j_temp_offset;
                                        decide_to_demand := True;                           {decide to demand since value>0}
                                    end; {if value>max value}
                            end; {j_temp <> i}
                    end; {j_temp}
            end;{contiguity = false}
    end; {Make_Demand}
  { ---------------------------------------------------------------  }
    procedure Make_Response_Decision;
        var
            payment: real;                              {expected payment, used in B3, B4}
            targets_upper_cost: real;               {cost to j if j rich enough, used in B3, B4}
            targets_cost_if_fight: real;                {expected damage to j, used in B3, B4}
    begin
        Calc_alliance(i, j, offset_target);     {determine W_att,W_def, aid_att, help_def. Also used in Conduct_Fight}
        case chrom[j, 2] of                                     {j's B rule for response}
            1: 
                agree_to_pay := false;                              {B1 never pay}
            2: 
                agree_to_pay := true;                               {B2 always pay}
            3:                                                          {B3 pay if cheaper for self than fighting}
                begin
                    payment := min(standard_demand, W[j]);          {expected payment}
                    if w_def = 0 then
                        targets_upper_cost := 0         {to avoid div by 0}
                    else
                        targets_upper_cost := destructiveness * W_att * (W[j] / W_def);
                    targets_cost_if_fight := min(targets_upper_cost, W[j]);
                    if targets_cost_if_fight > payment then  {considers only damage to j, not helpers}
                        Agree_To_Pay := True                    {pay if cheaper than fighting - for self (j)}
                    else
                        Agree_To_Pay := False;
{    if (year = 39) and (i = 5) and (j = 4) then}
{    begin}
{    writeln(' B3. payment=', payment : 5, 'targets_upper_cost= ', targets_upper_cost : 6 : 2);}
{    writeln('    targets_cost_if_fight= ', targets_cost_if_fight : 6 : 2, ' Wj=', W[j] : 6 : 2, ' Wdef=', W_def : 6 : 2);}
{end;}
{XXX }
                end;{case B3}
            4:                                                          {B4 fight if cheaper for self than paying}
                begin
                    payment := min(standard_demand, W[j]);          {expected payment}
                    if w_def = 0 then
                        targets_upper_cost := 0         {to avoid div by 0}
                    else
                        targets_upper_cost := destructiveness * W_att * (W[j] / W_def);
                    targets_cost_if_fight := min(targets_upper_cost, W[j]);
                    if targets_cost_if_fight >= payment then  {considers only damage to j, not helpers}
                        Agree_To_Pay := True                    {fight if cheaper than paying - for self (j)}
                    else
                        Agree_To_Pay := False;
                end;{case B4}
            5:                                                          {B5 pay if  fight_cost > toughness*payment}
                begin
                    payment := min(standard_demand, W[j]);          {expected payment}
                    if w_def = 0 then
                        targets_upper_cost := 0         {to avoid div by 0}
                    else
                        targets_upper_cost := destructiveness * W_att * (W[j] / W_def);
                    targets_cost_if_fight := min(targets_upper_cost, W[j]);
                    if targets_cost_if_fight > toughness * payment then  {considers only damage to j, not helpers}
                        Agree_To_Pay := True                    {pay if  fight_cost > toughness*payment - for self (j)}
                    else
                        Agree_To_Pay := False;
                end;{case B5}
            6:                                                          {B6 pay or fight at random}
                begin
                    if random < 0 then
                        Agree_To_Pay := True
                    else
                        Agree_To_Pay := False;
                end;{case B6}
        end;{case of chrom, ie j's B rule}
    end;{procedure Make_Response_Decision}
  { ---------------------------------------------------------------  }
    procedure Increase_Loyalty (from_person, to_person: actor_type);
        var
            temp_loyalty_pc: integer;           {temporary loyalty in percentage}
            loyalty_to_be_increased: boolean;       {who may increase loyalty: 0=all, 1=rare, 2=common}
{       dont need to worry about decreasing since it's say 0 if never increased}
    begin
        loyalty_to_be_increased := False;  {assume conditions won't be met}
        case loyalty_option of
            0:                              {everyone increases }
                loyalty_to_be_increased := True;
            1:                               {increase only if either from or to person is rare}
                if rare[from_person] or rare[to_person] then
                    loyalty_to_be_increased := True;
            2:                              {increase only if neither from or to person is rare}
                if not rare[from_person] and not rare[to_person] then
                    loyalty_to_be_increased := True;
        end; {case}
        if loyalty_to_be_increased then     {increase in normal way}
            begin
                temp_loyalty_pc := loyalty_pc[from_person, to_person] + loyalty_pc_increment;
                if temp_loyalty_pc > 100 then
                    Loyalty_pc[from_person, to_person] := 100           {max allowable}
                else
                    Loyalty_pc[from_person, to_person] := temp_loyalty_pc;
{test XXX}
{    writeln('IncLoy: newL, from,to ', loyalty_pc[from_person, to_person] : 7, from_person, to_person);}
            end; {loyalty option}
    end;{Increase_Loyalty}
 { ---------------------------------------------------------------  }
    procedure Decrease_Loyalty (from_person, to_person: actor_type);
        var
            temp_loyalty_pc: integer;           {temporary loyalty in percentage}
    begin
        temp_loyalty_pc := loyalty_pc[from_person, to_person] - loyalty_pc_increment;
        if temp_loyalty_pc < 0 then
            loyalty_pc[from_person, to_person] := 0         {min allowable}
        else
            Loyalty_pc[from_person, to_person] := temp_loyalty_pc;
{test XXX}
{writeln('DecLoy: newL, from,to ', loyalty_pc[from_person, to_person] : 7, from_person, to_person);}
    end;{Decrease_Loyalty}
 { ---------------------------------------------------------------  }
    procedure Make_Payment;
    begin
        if W[j] < Standard_Demand then                          {payment doesn't exceed wealth of payor}
            Tribute := W[j]
        else
            Tribute := Standard_Demand;
        W[i] := W[i] + Tribute;                             {j automatically pays}
        W[j] := W[j] - Tribute;
{writeln('XXX increase loyalty due to tribute', i, j);}
        Increase_Loyalty(i, j);                                 {receiver of tribute becomes more loyal to payer}
        Increase_Loyalty(j, i);                                 {payer of tribute becomes more loyal to receiver}
        trib_count[j, i] := trib_count[j, i] + 1;   {update total count from row to col, ie j to i}
        trib_periodic_count[j, i] := trib_periodic_count[j, i] + 1; {update per. count from row to col, ie j to i}
        if (report_person = i) or (report_person = j) then
            begin
                Agree_to_Pay := true;                   {tribute paid}
                Make_Person_Report;
            end;{if}
{write('XXX Pay: Year=', year : 3, ' to  i=', i : 3, ' from j=', j : 3);}
{writeln('  W[i]=', W[i] : 7 : 1, '     W[j]=', W[j] : 7 : 1);}
    end;
  { ---------------------------------------------------------------  }
    procedure Conduct_FIght;
        type
            report_string_type = packed array[1..imax] of char;
        var
            Damage_by_Attacker: real;                   {damage done by attacking alliance}
            Damage_by_Defender: real;
            k, L: actor_type;                               {generic actors}
            rebellion: boolean;                             {whether there is a potential rebellion}
            report_string: report_string_type;      {to report events with potential rebellion}
    begin
{writeln('Fight: Attacker= ', i : 2, '  Defender= ', j : 2);}
{writeln('  Wealth'' s of a, d before = ', W[i] : 6 : 2, '  ', W[j] : 6 : 2);}
        rebellion := false;                             {assume its not a potential rebellion}
        if report_rebellion and (loyalty_pc[i, j] = 100) and (w[j] > w[i]) then
            begin {xxx-         report a potential rebellion}
                rebellion := true;                          {to get a report after the fight}
                writeln('XXX fight: pop= ', pop : 3, ' year = ', year : 4, '. fight i , j ', i : 5, ' ', j : 5, ' %loy ( i , j ) = ', loyalty_pc[i, j] : 6);
                writeln('  attack of bigger target: Wi= ', W[i] : 6 : 1, ' Wj= ', W[j] : 6 : 1, ' W_def= ', W_def : 6 : 1, 'W_att = ', W_att : 6 : 1);
                write('W before fight:', w[1] : 6 : 1, '    ', w[2] : 6 : 1, '  ', w[3] : 6 : 1, '  ', w[4] : 6 : 1, '  ', w[5] : 6 : 1, '  ', w[6] : 6 : 1, '  ');
                writeln(w[7] : 6 : 1, ' ', w[8] : 6 : 1, '  ', w[9] : 6 : 1, '  ', w[10] : 6 : 1);
                for k := 1 to imax do                   {do row as outer loop}
                    begin
                        write(pop : 4, '    ', year : 4, '  ', k : 4, ' ');
                        for L := 1 to imax do
                            begin
                                write(loyalty_pc[k, L] : 4, '   ');
                                if L = imax then
                                    writeln;
                            end;{L}
                    end;{k}
            end;{report rebellion}
        Damage_by_Attacker := min(destructiveness * W_att, W_def);          {Can't do more than W_def}
        Damage_by_Defender := min(destructiveness * W_def, W_att);      {Can't do more than W_att}
        for k := 1 to imax do               {calc everyone's new wealth and loyalty}
            begin
                if aid_att[k] then                  {k is part of attack}
                    begin
                        if W_att > 0.0 then
                            W[k] := W[k] - Damage_by_Defender * loyalty_pc[k, i] * (W[k] / W_att) / 100;
                        for L := 1 to imax do
                            begin
                                if aid_att[L] then                                      {k, L on same side of fight}
                                    begin
                                        increase_loyalty(k, L);
                                    end;
                                if help_def[L] then                                 {k, L on different sides of fight}
                                    decrease_loyalty(k, L);
                            end;{L}
                    end;{k is in attack}
                if help_def[k] then             {k is part of defense}
                    begin
                        if W_def > 0.0 then
                            begin
                                W[k] := W[k] - Damage_by_Attacker * loyalty_pc[k, j] * (W[k] / W_def) / 100;
{write('damage_by_attacker, loyalty%, k, j, Wk new,W_def');}
 {XXX test}
{writeln(Damage_by_Attacker : 6 : 1, loyalty_pc[k, j] : 6, k : 3, j : 3, W[k] : 6 : 1, W_def : 6 : 1);  }
{XXX test}
                            end;
                        for L := 1 to imax do
                            begin
                                if help_def[L] then                                 {k, L on same side of fight}
                                    begin
{writeln('XXX from conduct fight: on same side,both def:', k, l);}
                                        increase_loyalty(k, L);
                                    end;
                                if aid_att[L] then                                      {k, L on different sides of fight}
                                    decrease_loyalty(k, L);
                            end;{L}
                    end;{K is in defense}
            end;{k}
        if (aid_att[report_person] or help_def[report_person]) then
            make_person_report;
        fights_since_last_report := fights_since_last_report + 1;       {count number of fights since periodic report}
        if rebellion then           {report string}
            begin
                for k := 1 to imax do
                    begin
                        report_string[k] := '-';                                {set to . if nothing applies}
                        if aid_att[k] then
                            report_string[k] := 'a';                                {a for aid attacker}
                        if help_def[k] then
                            report_string[k] := 'd';                                {d for help defender}
                        if k = i then
                            report_string[k] := 'A';                            {A for attacker}
                        if k = j then
                            report_string[k] := 'D';                            {D for defender}
                    end;{k}
                writeln(' Previous event= ', report_string : 12);
            end;{rebellion}
    end;{conduct fight}
  { ---------------------------------------------------------------  }
    procedure Produce_Wealth;
        var
            i: integer;
    begin
        for i := 1 to imax do
            begin
                Income := Productivity[i];
                W[i] := W[i] + Income;
                W_if_isolated[i] := W_if_isolated[i] + Income;
{Writeln('Year=', year : 3, '  i=', i : 3, '   Income=', Income : 7 : 1, ' Wealth=', W[i] : 7 : 1);}
            end;{production cycle}
    end;
  { ---------------------------------------------------------------  }
    procedure Report_Final_Output;
        var
            i, j: actor_type;
    begin
        for i := 1 to imax do
            begin
                writeln(final_wealth, pop : 3, '    ', rare[i] : 3, '   ', i : 3, ' ', W[i] : 9 : 1, '  ', W_if_isolated[i] : 9 : 1);
            end;
        writeln(final_wealth, 'Loyalty, in percent. From row to col.');
        writeln(final_wealth, 'pop  row Final W 1   2   3   4   5   6   7   8   9   10');
        for i := 1 to imax do                   {do row as outer loop}
            begin
                write(final_wealth, pop : 4, '  ', i : 4, ' ', W[i] : 9 : 1, '  ');
                for j := 1 to imax do
                    begin
                        write(final_wealth, loyalty_pc[i, j] : 4, ' ');
                        if j = imax then
                            writeln(final_wealth);
                    end;{j}
            end;{i}
        writeln(final_wealth, 'Count of  times trib paid by row to col.');
        writeln(final_wealth, 'pop  row Final W 1   2   3   4   5   6   7   8   9   10');
        for i := 1 to imax do                   {do row as outer loop}
            begin
                write(final_wealth, pop : 4, '  ', i : 4, ' ', W[i] : 9 : 1, '  ');
                for j := 1 to imax do
                    begin
                        write(final_wealth, trib_count[i, j] : 4, ' ');
                        if j = imax then
                            writeln(final_wealth);
                    end;{j}
            end;{i}
    end;
  { ---------------------------------------------------------------  }
    procedure Conduct_and_Report_Adaptation;                {This is "selection" procedure}
        var                                                                 {NOTE: mutation of rare rule is in initialize_pop}
            i_temp: actor_type;
            invade: integer;                                            {1 if invasion criterion met, 0 if not}
            gene_temp: gene_type;
            ave_common_wealth: real;
            ave_common_LR: real;                                {Left+Right range}
            rare_LR: integer;
            ave_common_successes: real;                         {number Wealth > Wealth if Iso}
            rare_success: integer;

    begin
        ave_common_wealth := 0;
        ave_common_LR := 0;
        ave_common_successes := 0;
        for i_temp := 1 to imax do
            begin
                if i_temp <> rare_indiv then            {cumulate stats for common indivs}
                    begin
                        ave_common_wealth := ave_common_wealth + w[i_temp];
{ave_common_LR := ave_common_LR + left_d[i_temp] + right_d[i_temp];}
{obsolete}
                        if W[i_temp] > W_if_isolated[i_temp] then
                            ave_common_successes := ave_common_successes + 1;
                    end; {if}
            end;{i_temp}
        ave_common_wealth := ave_common_wealth / (imax - 1);        {there are imax-1 common indivs}
        ave_common_LR := ave_common_LR / (imax - 1);
        ave_common_successes := ave_common_successes / (imax - 1);
{rare_LR := left_d[rare_indiv] + right_d[rare_indiv];}
{obsolete}
        if W[rare_indiv] > W_if_isolated[rare_indiv] then
            rare_success := 1
        else
            rare_success := 0;
        invade := 0;                    {assume unless following test works}
        duration_since_invasion := duration_since_invasion + 1;         {tally since last invasion}
        if W[rare_indiv] > ave_common_wealth then
            begin
                invade := 1;                                                {invade criterion; write line of invaders}
                invaders_so_far := invaders_so_far + 1;         {tally of how many times there has been an invasion}
                write(invaders_wealth, invaders_so_far : 4, '   ', duration_since_invasion : 4, '   ', pop : 4);
                write(invaders_wealth, '    ', chrom[common_indiv, 1] : 3, '    ', chrom[common_indiv, 2] : 2);
                write(invaders_wealth, '    ', chrom[common_indiv, 3] : 2, '    ', chrom[common_indiv, 4] : 2);
                write(invaders_wealth, '    ', chrom[rare_indiv, 1] : 3, '  ', chrom[rare_indiv, 2] : 2);
                write(invaders_wealth, '    ', chrom[rare_indiv, 3] : 2, '  ', chrom[rare_indiv, 4] : 2);
                write(invaders_wealth, '    ', ave_common_wealth : 8 : 2, ' ', w[rare_indiv] : 8 : 2);
                write(invaders_wealth, '    ', ave_common_LR : 7 : 2, ' ', rare_LR : 4);
                writeln(invaders_wealth, '  ', ave_common_successes : 7 : 2, '  ', rare_success : 4);
                duration_since_invasion := 0;                       {restart count of pops since last invasion}
            end;
{                   write line of adapt file whether or not invade  }
        write(adapt_wealth, pop : 4, '  ', chrom[common_indiv, 1] : 3, '    ', chrom[common_indiv, 2] : 2);
        write(adapt_wealth, '   ', chrom[common_indiv, 3] : 2, '    ', chrom[common_indiv, 4] : 2);
        write(adapt_wealth, '   ', chrom[rare_indiv, 1] : 3, '  ', chrom[rare_indiv, 2] : 2);
        write(adapt_wealth, '   ', chrom[rare_indiv, 3] : 2, '  ', chrom[rare_indiv, 4] : 2);
        write(adapt_wealth, '   ', ave_common_wealth : 8 : 2, ' ', w[rare_indiv] : 8 : 2, ' ', invade : 2);
        write(adapt_wealth, '   ', ave_common_LR : 7 : 2, ' ', rare_LR : 4);
        writeln(adapt_wealth, ' ', ave_common_successes : 7 : 2, '  ', rare_success : 4);
        if invade = 1 then                                              {make new common rules the old rare rule}
            begin
                for i_temp := 1 to imax do
                    begin
                        if i_temp <> rare_indiv then
                            begin
                                for gene_temp := 1 to gene_max do
                                    begin
                                        chrom[i_temp, gene_temp] := chrom[rare_indiv, gene_temp];
                                    end;{gene_temp}
                            end;{if i_temp<>rare_indiv}
                    end;{i_temp}
            end;{if invade}
    end;{Conduct_and_Report_Adaptation}
  { ---------------------------------------------------------------  }
  { ---------------------------------------------------------------  }
{M A I N    P R O G R A M }
begin
    initialize_run;                                                             {initialize whole run}
    for pop := 1 to pop_max do                                      {POPULATION CYCLE}
        begin
            Initialize_pop;                                                         {INITIALIZE current population}
            for year := 1 to year_max do                                 {YEARLY CYCLE}
                begin
                    for demand_phase := 1 to demand_phase_max do        {DEMAND CYCLE}
                        begin
                            Select_Active_Actor;
                            Make_Demand;
                            if Decide_To_Demand then                                {demand was made of j}
                                begin
                                    Make_Response_Decision;
                                    if Agree_To_Pay then
                                        Make_Payment
                                    else
                                        Conduct_Fight;
                                end;{decide to demand}
                        end; {demand cycle}
                    Produce_Wealth;                                      {PRODUCTION}
                    years_since_last_report := years_since_last_report + 1;
                    if years_since_last_report = periodic_report_freq then  {time to give periodic_output line}
                        Report_Periodic_Output;
                end; {yearly cycle}
            Report_Final_Output;                            {of a pop for a line of final_output}
            if adapt then
                Conduct_and_Report_Adaptation;          {of a pop for a line of adapt_output}
        end;{pop cycle}
    gettime(end_datetime);
    end_hour := end_datetime.hour;          {to force long interger}
    end_time := 60 * 60 * end_hour + 60 * end_datetime.minute + end_datetime.second;
    duration := end_time - start_time;
    Writeln('Duration of this run is ', duration : 5, ' seconds.');
    Writeln(final_wealth, 'Duration of this run is ', duration : 5, ' seconds.');
end.{main program}