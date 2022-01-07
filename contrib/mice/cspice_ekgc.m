%-Abstract
%
%   CSPICE_EKGC returns an element of string (character) data from a
%   specified row in a specified column of the set of rows matching
%   the previous cspice_ekfind SELECT query.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY,
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING,
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      selidx   the index for a column of interest satisfying the SELECT
%               clause, the column indices range from 1 to number of
%               columns in the SELECT clause.
%
%               [1,1] = size(selidx); int32 = class(selidx)
%
%      row      the index for a row in the column identified by 'selidx',
%               the column indices range from 1 to 'nmrows' where 'nmrows'
%               equals the total number of rows satisfying the SELECT clause.
%
%               [1,1] = size(row); int32 = class(row)
%
%      elment   the index for an element of the data at the 'selidx','row'
%               position; a scalar value at 'selidx', 'row' has 'elment'
%               value one.
%
%               [1,1] = size(elment); int32 = class(elment)
%
%      cdatln   the maximum length of the `cdata' output string.
%
%               [1,1] = size(cdatln); int32 = class(cdatln)
%
%   the call:
%
%      [cdata, null, found] = cspice_ekgc( selidx, row, elment, cdatln )
%
%   returns:
%
%      cdata    the string value of the requested element at data
%               location 'selidx', 'row', 'elment'.
%
%               [1,c2] = size(cdata); char = class(cdata)
%
%      null     a boolean indicating if 'cdata' has a null value.
%
%               [1,1] = size(null); logical = class(null)
%
%      found    a boolean indicating whether the specified value at
%               'selidx', 'row', 'elment' was found.
%
%               [1,1] = size(found); logical = class(found)
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Perform a query on an EK file that contains a database with
%      the Supplementary Engineering Data Records of the Viking Project
%      in order to retrieve the IMAGE_ID values (character strings)
%      that correspond to the images with IMAGE_NUMBER smaller than
%      a given value, ordered by IMAGE_NUMBER.
%
%
%      Use the EK kernel below to load the information from the
%      original Supplementary Engineering Data Record (SEDR) data
%      set generated by the Viking Project.
%
%         vo_sedr.bdb
%
%
%      Example code begins here.
%
%
%      function ekgc_ex1()
%         %
%         % Assign an EK file to load and a max string size
%         % for the cspice_ekgc return string.
%         %
%         EK     = 'vo_sedr.bdb';
%         MAXSTR = 1025;
%
%         %
%         % Load the EK.
%         %
%         cspice_furnsh( EK )
%
%         %
%         % The table 'VIKING_SEDR_DATA' has a column 'IMAGE_ID'
%         % of scalar strings.
%         %
%         % Define a set of constraints to perform a query on all
%         % loaded EK files (the SELECT clause). In this case select
%         % the column 'IMAGE_ID' from table 'VIKING_SEDR_DATA'
%         % sorted by 'IMAGE_NUMBER'.
%         %
%         query = ['Select IMAGE_ID from VIKING_SEDR_DATA '                ...
%                  'where IMAGE_NUMBER < 25860000  order by IMAGE_NUMBER'];
%
%         %
%         % Query the EK system for data rows matching the
%         % SELECT constraints.
%         %
%         [nmrows, error, errmsg] = cspice_ekfind( query );
%
%         %
%         % Check whether an error occurred while processing the
%         % SELECT clause. If so, output the error message.
%         %
%         if ( error )
%            printf( 'SELECT clause error: %s\n', errmsg );
%         end
%
%         %
%         % Loop over each row found matching the query.
%         %
%         for rowno = 1:nmrows
%
%            %
%            % Fetch the character data. We know the query returned
%            % one column and the column contains only scalar data,
%            % so the index of all elements is 1.
%            %
%            selidx = 1;
%            eltidx = 1;
%
%            %
%            % Use cspice_ekgc to retrieve the string from
%            % row/column position.
%            %
%            [cdata, isnull, found] = cspice_ekgc( selidx, rowno,          ...
%                                                  eltidx, MAXSTR );
%
%            %
%            % Output the value, if non-null data exist at the
%            % requested position.
%            %
%            if  ~isnull
%               fprintf( 'Row %3d: Character data: %s\n', rowno, cdata );
%            end
%
%         end
%
%         %
%         % Clear the kernel pool and database. Note, you don't normally
%         % unload an EK after a query, rather at the end of a program.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Row   1: Character data: 168C09
%      Row   2: Character data: 168C10
%      Row   3: Character data: 168C11
%      Row   4: Character data: 168C12
%      Row   5: Character data: 169C01
%      Row   6: Character data: 169C02
%      Row   7: Character data: 169C03
%      Row   8: Character data: 169C04
%      Row   9: Character data: 169C05
%      Row  10: Character data: 169C09
%      Row  11: Character data: 169C11
%      Row  12: Character data: 169C19
%      Row  13: Character data: 169C23
%      Row  14: Character data: 169C25
%      Row  15: Character data: 169C26
%      Row  16: Character data: 169C30
%      Row  17: Character data: 169C32
%      Row  18: Character data: 169C33
%      Row  19: Character data: 169C37
%      Row  20: Character data: 169C39
%      Row  21: Character data: 169C40
%      Row  22: Character data: 169C44
%      Row  23: Character data: 169C46
%      Row  24: Character data: 169C47
%      Row  25: Character data: 169C51
%      Row  26: Character data: 169C53
%
%
%-Particulars
%
%   Suppose a SELECT clause return data consisting of three columns (N=3)
%   and four rows (M=4):
%
%              col 1    col 2    col 3
%
%      row 1   val_11   val_12   val_13
%      row 2   val_21   val_22   val_23
%      row 3   val_31   val_32   val_33
%      row 4   val_41   val_42   val_43
%
%   with "col 2" and "col 3" containing scalar string data and
%   "val_42" containing a vector of K strings.
%
%   Retrieving the data elements depends on the values for the index set
%   "selidx," "row," and "elment."
%
%   Use the set
%
%      'selidx' = 2, 'row' = 3, 'elment' = 1
%
%   to fetch scalar "val_32."
%
%   Use the set
%
%      'selidx' = 3, 'row' = 4, 'elment' = 1
%
%   to fetch scalar "val_43."
%
%   Use the set
%
%      'selidx' = 2, 'row' = 4, 'elment' = K
%
%   to fetch the final element of vector "val_42"
%
%   `elment' is allowed to exceed the number of elements in the column
%   entry; if it does, `found' returns as false. This allows the caller
%   to read data from the column entry in a loop without checking the
%   number of available elements first.
%
%-Exceptions
%
%   1)  If the input argument `elment' is less than 1, the error
%       SPICE(INVALIDINDEX) is signaled by a routine in the call tree
%       of this routine and `found' is returned false. However, `elment'
%       is allowed to be greater than the number of elements in the
%       specified column entry; this allows the caller to read data
%       from the column entry in a loop without checking the number of
%       available elements first. If `elment' is greater than the number
%       of available elements, `found' is returned false.
%
%   2)  If `selidx' is outside of the range established by the last query
%       passed to the EK search engine, the error SPICE(INVALIDINDEX) is
%       signaled by a routine in the call tree of this routine and `found' is
%       returned false.
%
%   3)  If the input argument `row' is less than 1 or greater than the
%       number of rows matching the query, the error
%       SPICE(INVALIDINDEX) is signaled by a routine in the call tree
%       of this routine and `found' is returned false.
%
%   4)  If the specified column does not have character type, the
%       error SPICE(INVALIDTYPE) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If this routine is called when no E-kernels have been loaded,
%       the error SPICE(NOLOADEDFILES) is signaled by a routine in the
%       call tree of this routine.
%
%   6)  If any of the input arguments, `selidx', `row', `elment' or
%       `cdatln', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   7)  If any of the input arguments, `selidx', `row', `elment' or
%       `cdatln', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   This routine reads binary "sequence component" EK files.
%   In order for a binary EK file to be accessible to this routine,
%   the file must be "loaded" via a call to the routine cspice_furnsh.
%
%   Text format EK files cannot be used by this routine; they must
%   first be converted by binary format by the NAIF Toolkit utility
%   SPACIT.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   EK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.3.0, 10-AUG-2021 (EDW) (JDR)
%
%       Changed the input argument name "lenout" to "cdatln" for consistency
%       with other routines.
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement and and example's EK. Updated example
%       code to work with provided EK.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.2.1, 03-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.2.0, 10-MAY-2011 (EDW)
%
%       "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.0, 10-APR-2010 (EDW)
%
%-Index_Entries
%
%   fetch element from character column entry
%
%-&

function [cdata, null, found] = cspice_ekgc( selidx, row, elment, cdatln )

   switch nargin
      case 4

         selidx = zzmice_int(selidx);
         row    = zzmice_int(row);
         elment = zzmice_int(elment);
         cdatln = zzmice_int(cdatln);

      otherwise

         error ( [ 'Usage: [ `cdata`, null, found] = ' ...
                   'cspice_ekgc( selidx, row, elment, cdatln )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [cdata, null, found] = mice('ekgc_c', selidx, row, elment, cdatln );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      null  = zzmice_logical(null);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end


