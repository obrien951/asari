#ifndef _ASARI_BASE64_
#define _ASARI_BASE64_
#include <vector>

namespace asaristc {

class b64_decoder {
public:
  b64_decoder();

  //  void decode( unsigned char * &inp, unsigned char * &outp);
  void decode_base64(char *bin_data, const char *base64_data,
                     long unsigned int &bytes_to_decode);

protected:
  // Checks to see if the processor is bigendian or littlendian then forms
  //   the appropriate look-up table. The table isn't hard-coded to lower the
  //   memory footprint
  void form_table();

  bool tables_initiated_ = false;
  bool is_little_endian_;

  /* Look up tables for base64 to binary conversion */
  /* Copy values out of this table 3 bits at a time.
   * Decodes the bulk of the data */
  std::vector<unsigned char> main_lookup_;

  /* Dealing with padded data. Are used to decode the first second and third
   * dangling characters respectively */
  std::vector<unsigned char> back_lookup_1_;
  std::vector<unsigned char> back_lookup_2_;
  std::vector<unsigned char> back_lookup_3_;

  unsigned int base_table_[128] = {
      // basic base64 charset table
      0,  //  NUL
      0,  //  SOH
      0,  //  STX
      0,  //  ETX
      0,  //  EOT
      0,  //  ENQ
      0,  //  ACK
      0,  //  BEL
      0,  //   BS
      0,  //   HT
      0,  //   LF
      0,  //   VT
      0,  //   FF
      0,  //   CR
      0,  //   SO
      0,  //   SI
      0,  //  DLE
      0,  //  DC1
      0,  //  DC2
      0,  //  DC3
      0,  //  DC4
      0,  //  NAK
      0,  //  SYN
      0,  //  ETB
      0,  //  CAN
      0,  //   EM
      0,  //  SUB
      0,  //  ESC
      0,  //   FS
      0,  //   GS
      0,  //   RS
      0,  //   US
      0,  //   SP
      0,  //    !
      0,  //    "
      0,  //    #
      0,  //    $
      0,  //    %
      0,  //    &
      0,  //    '
      0,  //    (
      0,  //    )
      0,  //    *
      62, //    +
      0,  //    ,
      0,  //    -
      0,  //    .
      63, //    /
      52, //    0,
      53, //    1
      54, //    2
      55, //    3
      56, //    4
      57, //    5
      58, //    6
      59, //    7
      60, //    8
      61, //    9
      0,  //    :
      0,  //    ;
      0,  //    <
      0,  //    =
      0,  //    >
      0,  //    ?
      0,  //    @
      0,  //    A
      1,  //    B
      2,  //    C
      3,  //    D
      4,  //    E
      5,  //    F
      6,  //    G
      7,  //    H
      8,  //    I
      9,  //    J
      10, //    K
      11, //    L
      12, //    M
      13, //    N
      14, //    O
      15, //    P
      16, //    Q
      17, //    R
      18, //    S
      19, //    T
      20, //    U
      21, //    V
      22, //    W
      23, //    X
      24, //    Y
      25, //    Z
      0,  //    [
      0,  //    '\'
      0,  //    ]
      0,  //    ^
      0,  //    _
      0,  //    `
      26, //    a
      27, //    b
      28, //    c
      29, //    d
      30, //    e
      31, //    f
      32, //    g
      33, //    h
      34, //    i
      35, //    j
      36, //    k
      37, //    l
      38, //    m
      39, //    n
      40, //    o
      41, //    p
      42, //    q
      43, //    r
      44, //    s
      45, //    t
      46, //    u
      47, //    v
      48, //    w
      49, //    x
      50, //    y
      51, //    z
      0,  //    {
      0,  //    |
      0,  //    }
      0,  //    ~
      0   //  DEL
  };
};

} // namespace asaristc

#else
#endif
