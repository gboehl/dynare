% Test translation from matched_moments block to M_.matched_moments

@#define estimParams = 0
@#define MoM_Method = "SMM"
@#define NoEstim = 1
@#define Extended_Matched_Moments_Checks
@#include "AnScho_MoM_common.inc"

% get indices in declaration order
iYGR  = strmatch('YGR',  M_.endo_names,'exact');
iINFL = strmatch('INFL', M_.endo_names,'exact');
iINT  = strmatch('INT',  M_.endo_names,'exact');

% M_.matched_moments has the following structure:
%   - first entry: index number of variable in declaration order
%   - second entry: lead or lag
%   - third entry: power

matched_moments_orig = {

    %first-order product moments
    [iYGR       ]  [0 ],  [1];
    [iYGR       ]  [-1],  [1];
    [iYGR       ]  [0 ],  [1];
    [iYGR       ]  [0 ],  [1];
    [iYGR       ]  [0 ],  [1];
    [iINFL      ]  [2 ],  [1];
    [iINT       ]  [-2],  [1];
    [iINT       ]  [2 ],  [1];

    %second-order contemporenous product moments
    [iYGR       ]  [0     ],  [2  ];
    [iYGR       ]  [-1    ],  [2  ];
    [iYGR  iYGR ]  [0    0],  [1 1];
    [iYGR  iYGR ]  [1    1],  [1 1];
    [iYGR  iINFL]  [0    0],  [1 1];
    [iYGR  iINFL]  [0    0],  [1 1];
    [iINT  iYGR ]  [2    2],  [1 1];
    [iINT  iYGR ]  [2    2],  [1 1];
    [iINT  iYGR ]  [2    2],  [1 1];
    [iINT  iYGR ]  [2    2],  [1 1];
    [iINFL iINFL]  [-2  -2],  [1 1];
    [iINFL iINT ]  [-1  -1],  [1 1];
    [iINT  iINT ]  [0    0],  [1 1];
    [iINT  iINT ]  [0    0],  [1 1];

    %second-order temporal product moments
    [iYGR  iYGR ]  [0  -1],  [1 1];
    [iYGR  iYGR ]  [2   3],  [1 1];
    [iYGR  iYGR ]  [2   3],  [1 1];
    [iYGR  iYGR ]  [-1 -2],  [1 1];
    [iINT  iINT ]  [-2  3],  [1 1];
    [iINFL iINFL]  [0   2],  [1 1];
    [iYGR  iINFL]  [5  -3],  [1 1];
    [iYGR  iINFL]  [5  -3],  [1 1];
    [iINT  iINFL]  [2   3],  [1 1];

    [iYGR            ]  [0      ],  [3     ];
    [iYGR  iYGR      ]  [-3  -3 ],  [1  2  ];
    [iYGR  iYGR      ]  [-3  -3 ],  [1  2  ];
    [iYGR  iYGR  iYGR]  [0  0  0],  [1  1  1];
    [iINT  iINFL     ]  [-2  -1 ],  [2  4  ];
    [iINFL iINT      ]  [1   0 ],   [4  2  ];
    [iINFL iINT      ]  [1   0 ],   [4  2  ];
%    [iYGR  iINFL iINT]  [0  -3  5],   [2  4  6];
%    [iINFL iYGR  iINT]  [-3  0  5],   [4  2  6];
%    [iINFL iYGR  iINT]  [-3  0  5],   [4  2  6];
};

% Removed duplicate moment conditions
matched_moments_no_duplicate= {

    %first-order product moments
    [iYGR       ]  [0 ],  [1];
%    [iYGR       ]  [-1],  [1];
%    [iYGR       ]  [0 ],  [1];
%    [iYGR       ]  [0 ],  [1];
%    [iYGR       ]  [0 ],  [1];
    [iINFL      ]  [2 ],  [1];
    [iINT       ]  [-2],  [1];
%    [iINT       ]  [2 ],  [1];

    %second-order contemporenous product moments
    [iYGR       ]  [0     ],  [2  ];
%    [iYGR       ]  [-1    ],  [2  ];
%    [iYGR  iYGR ]  [0    0],  [1 1];
%    [iYGR  iYGR ]  [1    1],  [1 1];
    [iYGR  iINFL]  [0    0],  [1 1];
%    [iYGR  iINFL]  [0    0],  [1 1];
    [iINT  iYGR ]  [2    2],  [1 1];
%    [iINT  iYGR ]  [2    2],  [1 1];
%    [iINT  iYGR ]  [2    2],  [1 1];
%    [iINT  iYGR ]  [2    2],  [1 1];
    [iINFL iINFL]  [-2  -2],  [1 1];
    [iINFL iINT ]  [-1  -1],  [1 1];
    [iINT  iINT ]  [0    0],  [1 1];
%    [iINT  iINT ]  [0    0],  [1 1];

    %second-order temporal product moments
    [iYGR  iYGR ]  [0  -1],  [1 1];
%    [iYGR  iYGR ]  [2   3],  [1 1];
%    [iYGR  iYGR ]  [2   3],  [1 1];
%    [iYGR  iYGR ]  [-1 -2],  [1 1];
    [iINT  iINT ]  [-2  3],  [1 1];
    [iINFL iINFL]  [0   2],  [1 1];
    [iYGR  iINFL]  [5  -3],  [1 1];
%    [iYGR  iINFL]  [5  -3],  [1 1];
    [iINT  iINFL]  [2   3],  [1 1];

    [iYGR            ]  [0      ],  [3     ];
%    [iYGR  iYGR      ]  [-3  -3 ],  [1  2  ];
%    [iYGR  iYGR      ]  [-3  -3 ],  [1  2  ];
%    [iYGR  iYGR  iYGR]  [0  0  0],  [1  1  1];
    [iINT  iINFL     ]  [-2  -1 ],  [2  4  ];
%    [iINFL iINT      ]  [1   0 ],   [4  2  ];
%    [iINFL iINT      ]  [1   0 ],   [4  2  ];
%    [iYGR  iINFL iINT]  [0  -3  5],   [2  4  6];
%    [iINFL iYGR  iINT]  [-3  0  5],   [4  2  6];
%    [iINFL iYGR  iINT]  [-3  0  5],   [4  2  6];
};



if ~isequal(M_.matched_moments_orig,matched_moments_orig)
    error('Translation to matched_moments-block failed!')
else
    fprintf('Translation to matched_moments-block successful!\n\n')
end

if ~isequal(M_.matched_moments,matched_moments_no_duplicate)
    error('Removal of duplicate moment conditions failed!')
else
    fprintf('Removal of duplicate moment conditions was successful!\n\n')
end
