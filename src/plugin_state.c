/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <link.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

struct state_region {
    unsigned char *address;
    unsigned char *snapshot;
    size_t size;
};

static struct state_region *regions;
static size_t region_count;

static int locate_regions(struct dl_phdr_info *info, size_t size, void *data)
{
    const uintptr_t base = *(const uintptr_t *)data;
    const uintptr_t page = (uintptr_t)sysconf(_SC_PAGESIZE);
    uintptr_t relro_start = 0;
    uintptr_t relro_end = 0;
    (void)size;
    if ((uintptr_t)info->dlpi_addr != base) {
        return 0;
    }
    for (size_t index = 0; index < info->dlpi_phnum; ++index) {
        const ElfW(Phdr) *header = &info->dlpi_phdr[index];
        if (header->p_type == PT_GNU_RELRO) {
            relro_start =
                (info->dlpi_addr + header->p_vaddr) & ~(page - 1);
            relro_end =
                (info->dlpi_addr + header->p_vaddr + header->p_memsz +
                 page - 1) &
                ~(page - 1);
        }
    }
    for (size_t index = 0; index < info->dlpi_phnum; ++index) {
        const ElfW(Phdr) *header = &info->dlpi_phdr[index];
        if (header->p_type == PT_LOAD && (header->p_flags & PF_W) != 0) {
            const uintptr_t start = info->dlpi_addr + header->p_vaddr;
            const uintptr_t end = start + header->p_memsz;
            if (relro_start >= end || relro_end <= start) {
                ++region_count;
            } else {
                region_count += start < relro_start ? 1 : 0;
                region_count += relro_end < end ? 1 : 0;
            }
        }
    }
    regions = calloc(region_count, sizeof(*regions));
    if (regions == NULL) {
        region_count = 0;
        return 1;
    }
    size_t target = 0;
    for (size_t index = 0; index < info->dlpi_phnum; ++index) {
        const ElfW(Phdr) *header = &info->dlpi_phdr[index];
        if (header->p_type != PT_LOAD || (header->p_flags & PF_W) == 0) {
            continue;
        }
        const uintptr_t start = info->dlpi_addr + header->p_vaddr;
        const uintptr_t end = start + header->p_memsz;
        const bool overlaps = relro_start < end && relro_end > start;
        const uintptr_t starts[2] = {start, relro_end};
        const uintptr_t ends[2] = {overlaps ? relro_start : end, end};
        for (size_t part = 0; part < 2; ++part) {
            if (starts[part] >= ends[part] ||
                (part == 1 && !overlaps)) {
                continue;
            }
            regions[target].address = (unsigned char *)starts[part];
            regions[target].size = (size_t)(ends[part] - starts[part]);
            regions[target].snapshot = malloc(regions[target].size);
            if (regions[target].snapshot == NULL) {
                return 1;
            }
            ++target;
        }
    }
    return 1;
}

void potlibc_state_checkpoint(void)
{
    Dl_info image;
    uintptr_t base;
    if (dladdr((void *)&potlibc_state_checkpoint, &image) == 0) {
        abort();
    }
    base = (uintptr_t)image.dli_fbase;
    dl_iterate_phdr(locate_regions, &base);
    if (regions == NULL || region_count == 0) {
        abort();
    }
    for (size_t index = 0; index < region_count; ++index) {
        if (regions[index].snapshot == NULL) {
            abort();
        }
        memcpy(regions[index].snapshot, regions[index].address,
               regions[index].size);
    }
}

void potlibc_state_restore(void)
{
    struct state_region *saved_regions = regions;
    const size_t saved_count = region_count;
    for (size_t index = 0; index < saved_count; ++index) {
        memcpy(saved_regions[index].address, saved_regions[index].snapshot,
               saved_regions[index].size);
    }
}

void potlibc_state_release(void)
{
    for (size_t index = 0; index < region_count; ++index) {
        free(regions[index].snapshot);
    }
    free(regions);
}

static int host_abi_accepted;

void potlibc_plugin_note_host_abi(void)
{
    host_abi_accepted = 1;
}

int potlibc_plugin_host_abi_ok(void)
{
    return host_abi_accepted;
}
