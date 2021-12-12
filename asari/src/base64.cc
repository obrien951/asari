#include <cstring>
#include <vector>

#include "base64.h"

namespace asaristc {

b64_decoder::b64_decoder() { form_table(); }

void b64_decoder::form_table() {
  if (tables_initiated_) {
    return;
  }

  int i, j, k;

  int prob_int = 1;
  is_little_endian_ = (*((char *)&prob_int) > 0 ? true : false);

  back_lookup_1_.resize(32767);
  back_lookup_2_.resize(32767);
  back_lookup_3_.resize(32767);

  main_lookup_.resize(16777214);

  for (i = '+'; i <= 'z'; i++) {
    for (j = '+'; j <= 'z'; j++) {
      int dest_adr = (i << 8) | j; /* could be +, but this is commutative */
      back_lookup_1_[dest_adr] = (base_table_[i] << 2) | (base_table_[j] >> 4);
      back_lookup_2_[dest_adr] = (base_table_[i] << 4) | (base_table_[j] >> 2);
      back_lookup_3_[dest_adr] = (base_table_[i] << 6) | (base_table_[j]);
    }
  }

  if (is_little_endian_) {
    for (i = '+'; i <= 'z'; i++) {
      for (j = '+'; j <= 'z'; j++) {
        for (k = '+'; k <= 'z'; k++) {
          int mailman = 0;
          char *c = (char *)&mailman;
          *c++ = i;
          *c++ = j;
          *c = k;
          mailman *= 2;
          main_lookup_[mailman++] = back_lookup_1_[(i << 8) | j];
          main_lookup_[mailman++] = back_lookup_2_[(j << 8) | k];
        }
      }
    }
  }
}

void b64_decoder::decode_base64(char *bin_data, const char *base64_data,
                                long unsigned int &bytes_to_decode) {
  unsigned char *bin_dest = (unsigned char *)bin_data;
  const char *char_src = base64_data;
  int bytes_left = static_cast<int>(bytes_to_decode);
  unsigned short int f;
  int data_add = 0;
  while (bytes_left >= 3) {
    memcpy((char *)&data_add, char_src, 3);
    memcpy(bin_dest, &main_lookup_[data_add * 2], 2);
    *(bin_dest + 2) =
        back_lookup_3_[((*(char_src + 2)) << 8) | (*(char_src + 3))];
    char_src += 4;
    bin_dest += 3;
    bytes_left -= 3;
  }

  while (bytes_left--) {
    *bin_dest++ = back_lookup_1_[f = (((*char_src) << 8) | (*(char_src + 1)))];
    if (bytes_left--) {
      *bin_dest++ = back_lookup_2_[f = ((f << 8) | (*(char_src + 2)))];
      if (bytes_left--) {
        *bin_dest++ = back_lookup_3_[(short int)((f << 8) | (*(char_src + 3)))];
      } else {
        break;
      }
    } else {
      break;
    }
    char_src += 4;
  }
}

} // namespace asaristc
